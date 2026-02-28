/**
 * basic.test.ts — Integration test: real FM-Index → StrandView consumer.
 *
 * Builds an FM-Index from the synthetic test genome, streams alignments for a
 * query sequence known to be present, and validates the records via StrandView.
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { resolve } from 'path';
import {
  buildSchema,
  computeStrandMap,
  initStrandHeader,
  StrandView,
} from '@strand/core';

// The compiled native addon (requires: npm run build first).
// napi-rs exports NativeFmIndex (camelCase); NativeFMIndex is a TS type alias.
// eslint-disable-next-line @typescript-eslint/no-require-imports
const { NativeFmIndex } = require('../index.js') as {
  NativeFmIndex: new (path: string) => {
    streamAlignments(sab: Buffer, queries: string[], maxMismatches?: number): Promise<void>;
    chromNames(): string[];
  };
};
// Alias for readability in tests.
const NativeFMIndex = NativeFmIndex;

const GENOME_PATH = resolve(__dirname, 'fixtures/test_genome.fa');

// Alignment record schema — MUST match the Rust constants in strand.rs.
// buildSchema determines byte offsets by applying its own alignment rules:
//   chrom      utf8_ref → offset 0   (4-byte intern table handle)
//   pos        u32      → offset 4
//   query_id   u32      → offset 8
//   strand     bool8    → offset 12
//   mismatches u8       → offset 13
//   score      f32      → offset 16  (f32 requires 4-byte alignment; padded from 14)
//   record_stride = 20
const schema = buildSchema([
  { name: 'chrom',      type: 'utf8_ref' },
  { name: 'pos',        type: 'u32'      },
  { name: 'query_id',   type: 'u32'      },
  { name: 'strand',     type: 'bool8'    },
  { name: 'mismatches', type: 'u8'       },
  { name: 'score',      type: 'f32'      },
]);

const INDEX_CAPACITY = 65536;

function makeSab(chromNames: string[] = []) {
  const map = computeStrandMap({
    schema,
    index_capacity: INDEX_CAPACITY,
    heap_capacity: 0,
    query: { assembly: '', chrom: '', start: 0, end: 0 },
    estimated_records: 0,
  });
  const sab = new SharedArrayBuffer(map.total_bytes);
  initStrandHeader(sab, map, chromNames.length ? { intern: chromNames } : undefined);
  const view = new StrandView(sab);
  if (chromNames.length) view.updateInternTable(chromNames);
  return { sab, map, view };
}

// ─── Basic integration test ───────────────────────────────────────────────────

describe('NativeFMIndex basic integration', () => {
  let searcher: InstanceType<typeof NativeFMIndex>;

  beforeAll(() => {
    searcher = new NativeFMIndex(GENOME_PATH);
  });

  it('finds exact-match hits and writes them into the SAB', async () => {
    const chromNames = searcher.chromNames();
    const { sab, view } = makeSab(chromNames);

    // 'ATGATGATGATGATGATGATG' appears multiple times in test_genome.fa.
    const query = 'ATGATGATGATGATGATGATG';

    // Start streaming; await completion.
    await searcher.streamAlignments(Buffer.from(sab), [query]);

    // All records should now be committed.
    const committed = view.committedCount;
    expect(committed).toBeGreaterThan(0);

    const cursor = view.allocateCursor();
    for (let seq = 0; seq < committed; seq++) {
      const ok = cursor.seek(seq);
      expect(ok).toBe(true);

      // chrom is a utf8_ref handle — resolves to a chromosome name string.
      expect(typeof cursor.getRef('chrom')).toBe('string');
      expect(cursor.getU32('pos')).toBeGreaterThanOrEqual(0);
      expect(cursor.getU32('query_id')).toBe(0);
      expect(cursor.getBool('strand')).not.toBeNull();
      expect(cursor.getU8('mismatches')).toBe(0);   // exact match
      expect(cursor.getF32('score')).toBeCloseTo(1.0); // 1/(1+0) = 1.0
    }

    expect(view.status).toBe('eos');
  });

  it('returns zero records for a query not present in the genome', async () => {
    const { sab, view } = makeSab();

    await searcher.streamAlignments(Buffer.from(sab), ['NNNNNNNNNNNNNNNNNNNNNNNNNNNN']);

    expect(view.committedCount).toBe(0);
    expect(view.status).toBe('eos');
  });

  it('handles multiple queries in a single call', async () => {
    const { sab, view } = makeSab();

    await searcher.streamAlignments(Buffer.from(sab), [
      'ATGATGATGATGATGATGATG',
      'GCTAGCTAGCTAGCTAGCTAG',
    ]);

    expect(view.committedCount).toBeGreaterThan(0);
    expect(view.status).toBe('eos');
  });
});

// ─── Abort test ───────────────────────────────────────────────────────────────

describe('NativeFMIndex abort', () => {
  it('resolves cleanly when consumer signals abort mid-stream', async () => {
    // Use a tiny ring to force backpressure, making abort observable.
    const map = computeStrandMap({
      schema,
      index_capacity: 4,  // very small ring
      heap_capacity: 0,
      query: { assembly: '', chrom: '', start: 0, end: 0 },
      estimated_records: 0,
    });
    const sab = new SharedArrayBuffer(map.total_bytes);
    initStrandHeader(sab, map);

    const ctrl = new Int32Array(sab);
    const searcher = new NativeFMIndex(GENOME_PATH);

    // Set CTRL_ABORT = 1 (index 13) before the stream starts to ensure it's
    // seen immediately; a real mid-stream abort would set it after some records.
    Atomics.store(ctrl, 13, 1);

    // The stream should resolve (not hang or throw) with STATUS_ERROR.
    await expect(
      searcher.streamAlignments(Buffer.from(sab), ['ATGATGATGATGATGATGATG'])
    ).resolves.toBeUndefined();

    // Producer should have set STATUS_ERROR (= 3).
    const status = Atomics.load(ctrl, 10);
    expect(status).toBe(3); // STATUS_ERROR
  });
});

// ─── Backpressure test ────────────────────────────────────────────────────────

describe('NativeFMIndex backpressure', () => {
  it('stream completes with a tiny ring when consumer drains continuously', async () => {
    // Use a ring of capacity 4 — forces the producer to stall repeatedly.
    const map = computeStrandMap({
      schema,
      index_capacity: 4,
      heap_capacity: 0,
      query: { assembly: '', chrom: '', start: 0, end: 0 },
      estimated_records: 0,
    });
    const sab = new SharedArrayBuffer(map.total_bytes);
    initStrandHeader(sab, map);

    const ctrl = new Int32Array(sab);
    const view = new StrandView(sab);
    const searcher = new NativeFMIndex(GENOME_PATH);

    // Start streaming in the background.
    const streamPromise = searcher.streamAlignments(
      Buffer.from(sab),
      ['ATGATGATGATGATGATGATG'],
    );

    // Consumer drain loop: periodically advance READ_CURSOR so the producer
    // can make progress. Continue until STATUS_EOS or STATUS_ERROR is set.
    const drainLoop = async () => {
      const CTRL_COMMIT_SEQ_IDX = 8;
      const CTRL_READ_CURSOR_IDX = 9;
      const CTRL_STATUS_IDX = 10;

      let acked = 0;
      // Poll every 2ms, drain until terminal status.
      for (let i = 0; i < 2000; i++) {
        await new Promise(r => setTimeout(r, 2));
        const committed = Atomics.load(ctrl, CTRL_COMMIT_SEQ_IDX);
        if (committed > acked) {
          acked = committed;
          // Advance cursor to release ring slots.
          Atomics.store(ctrl, CTRL_READ_CURSOR_IDX, acked);
          Atomics.notify(ctrl, CTRL_READ_CURSOR_IDX, 1);
        }
        const status = Atomics.load(ctrl, CTRL_STATUS_IDX);
        if (status === 2 /* EOS */ || status === 3 /* ERROR */) break;
      }
    };

    // Run drain loop concurrently with the stream.
    await Promise.all([streamPromise, drainLoop()]);

    expect(['eos', 'error']).toContain(view.status);
  }, 15_000 /* generous timeout for slow CI */);
});
