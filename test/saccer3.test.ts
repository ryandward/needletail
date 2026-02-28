/**
 * saccer3.test.ts — Integration tests against the real SacCer3 genome.
 *
 * Sequences and expected hit counts are verified independently:
 *   - bowtie2 (exact end-to-end, 0 mismatches)
 *   - brute-force Python scan of sacCer3.fa
 *
 * TATTTATACCCATTCCCTCA — 1 hit, chrI (NC_001133.9) pos 2261 fwd
 * TGGGATTCCATTGTTGATAA — 80 hits, all 16 nuclear chroms, both strands
 *                         (Ty retrotransposon LTR-like repeat)
 */

import { describe, it, expect, beforeAll } from 'vitest';
import { resolve } from 'path';
import {
  buildSchema,
  computeStrandMap,
  initStrandHeader,
  StrandView,
} from '@strand/core';

// eslint-disable-next-line @typescript-eslint/no-require-imports
const { NativeFmIndex } = require('../index.js') as {
  NativeFmIndex: new (path: string) => {
    streamAlignments(sab: Buffer, queries: string[], maxMismatches?: number): Promise<void>;
    chromNames(): string[];
  };
};
const NativeFMIndex = NativeFmIndex;

const SACCER3 = resolve(__dirname, '../../../Git/SeqChain/tests/data/saccer3/sacCer3.fa');

const UNIQUE_QUERY   = 'TATTTATACCCATTCCCTCA'; // 1 exact hit in SacCer3
const REPEAT_QUERY   = 'TGGGATTCCATTGTTGATAA'; // 80 exact hits (Ty LTR repeat)

// ── Schema: chrom as utf8_ref (intern table handle), query_id for attribution ─

const schema = buildSchema([
  { name: 'chrom',      type: 'utf8_ref' },
  { name: 'pos',        type: 'u32'      },
  { name: 'query_id',   type: 'u32'      },
  { name: 'strand',     type: 'bool8'    },
  { name: 'mismatches', type: 'u8'       },
  { name: 'score',      type: 'f32'      },
]);

// makeSab populates the v4 metadata region with chromosome names so any
// consumer can resolve chrom handles without a side-channel call.
function makeSab(chromNames: string[]) {
  const map = computeStrandMap({
    schema,
    index_capacity: 65536,
    heap_capacity: 0,
    query: { assembly: '', chrom: '', start: 0, end: 0 },
    estimated_records: 0,
  });
  const sab = new SharedArrayBuffer(map.total_bytes);
  initStrandHeader(sab, map, { intern: chromNames });

  const view = new StrandView(sab);
  view.updateInternTable(chromNames);
  return { sab, map, view };
}

// ── Single shared index — building takes ~1.2 s for 12 MB genome ──────────────

describe('NativeFMIndex — SacCer3 genome', () => {
  let searcher: InstanceType<typeof NativeFMIndex>;
  let chromNames: string[];

  beforeAll(() => {
    searcher   = new NativeFMIndex(SACCER3);
    chromNames = searcher.chromNames();
  }, 30_000);

  // ── chromNames ───────────────────────────────────────────────────────────────

  it('chromNames() returns 17 NCBI accession names starting with NC_001133.9', () => {
    expect(chromNames).toHaveLength(17);
    expect(chromNames[0]).toBe('NC_001133.9');
    for (const name of chromNames) {
      expect(name).toMatch(/^NC_\d+\.\d+$/);
    }
  });

  // ── Single hit ──────────────────────────────────────────────────────────────

  it('finds exactly 1 hit for unique query TATTTATACCCATTCCCTCA', async () => {
    const { sab, view } = makeSab(chromNames);

    await searcher.streamAlignments(Buffer.from(sab), [UNIQUE_QUERY]);

    expect(view.committedCount).toBe(1);

    const cursor = view.allocateCursor();
    cursor.seek(0);

    // Verified against bowtie2: NC_001133.9 (chrom handle 0 = chrI), pos 2261, fwd.
    // cursor.getRef resolves the utf8_ref handle through the intern table.
    expect(cursor.getRef('chrom')).toBe('NC_001133.9');
    expect(cursor.getU32('pos')).toBe(2261);
    expect(cursor.getU32('query_id')).toBe(0);
    expect(cursor.getBool('strand')).toBe(true);
    expect(cursor.getU8('mismatches')).toBe(0);
    expect(cursor.getF32('score')).toBeCloseTo(1.0);

    expect(view.status).toBe('eos');
  });

  // ── Multiple hits ───────────────────────────────────────────────────────────

  it('finds exactly 80 hits for repeat query TGGGATTCCATTGTTGATAA', async () => {
    const { sab, view } = makeSab(chromNames);

    await searcher.streamAlignments(Buffer.from(sab), [REPEAT_QUERY]);

    expect(view.committedCount).toBe(80);

    const cursor = view.allocateCursor();
    for (let seq = 0; seq < 80; seq++) {
      cursor.seek(seq);
      expect(cursor.getU8('mismatches')).toBe(0);
      expect(cursor.getF32('score')).toBeCloseTo(1.0);
      expect(cursor.getU32('query_id')).toBe(0);
      // chrom handle resolves to an NCBI accession on one of the 16 nuclear chroms.
      const name = cursor.getRef('chrom');
      expect(name).toMatch(/^NC_\d+\.\d+$/);
    }

    expect(view.status).toBe('eos');
  });

  it('both strands represented in the 80-hit repeat result', async () => {
    const { sab, view } = makeSab(chromNames);

    await searcher.streamAlignments(Buffer.from(sab), [REPEAT_QUERY]);

    const cursor = view.allocateCursor();
    let fwd = 0, rev = 0;
    for (let seq = 0; seq < view.committedCount; seq++) {
      cursor.seek(seq);
      cursor.getBool('strand') ? fwd++ : rev++;
    }
    expect(fwd).toBeGreaterThan(0);
    expect(rev).toBeGreaterThan(0);
  });

  // ── Both queries in a single call ───────────────────────────────────────────

  it('returns 81 total hits and correct query_id attribution for both queries', async () => {
    const { sab, view } = makeSab(chromNames);

    await searcher.streamAlignments(Buffer.from(sab), [UNIQUE_QUERY, REPEAT_QUERY]);

    expect(view.committedCount).toBe(81);

    const cursor = view.allocateCursor();
    let fromUnique = 0, fromRepeat = 0;
    for (let seq = 0; seq < 81; seq++) {
      cursor.seek(seq);
      const qid = cursor.getU32('query_id');
      if (qid === 0) fromUnique++;
      else if (qid === 1) fromRepeat++;
    }
    expect(fromUnique).toBe(1);
    expect(fromRepeat).toBe(80);

    expect(view.status).toBe('eos');
  });

  // ── Mismatch search ──────────────────────────────────────────────────────────

  it('maxMismatches=1 finds the exact hit plus near-misses for unique query', async () => {
    const { sab, view } = makeSab(chromNames);

    await searcher.streamAlignments(Buffer.from(sab), [UNIQUE_QUERY], 1);

    expect(view.committedCount).toBeGreaterThanOrEqual(1);

    const cursor = view.allocateCursor();
    let exactHits = 0;
    for (let seq = 0; seq < view.committedCount; seq++) {
      cursor.seek(seq);
      const mm = cursor.getU8('mismatches');
      expect(mm).toBeLessThanOrEqual(1);
      expect(cursor.getF32('score')).toBeCloseTo(1.0 / (1.0 + mm));
      if (mm === 0) exactHits++;
    }
    expect(exactHits).toBeGreaterThanOrEqual(1);
    expect(view.status).toBe('eos');
  });

  it('maxMismatches=1 returns more hits than exact for repeat query', async () => {
    const { sab: sabExact, view: viewExact } = makeSab(chromNames);
    await searcher.streamAlignments(Buffer.from(sabExact), [REPEAT_QUERY], 0);

    const { sab: sabMm1, view: viewMm1 } = makeSab(chromNames);
    await searcher.streamAlignments(Buffer.from(sabMm1), [REPEAT_QUERY], 1);

    expect(viewMm1.committedCount).toBeGreaterThan(viewExact.committedCount);

    const cursor = viewMm1.allocateCursor();
    for (let seq = 0; seq < viewMm1.committedCount; seq++) {
      cursor.seek(seq);
      expect(cursor.getU8('mismatches')).toBeLessThanOrEqual(1);
    }
    expect(viewMm1.status).toBe('eos');
  });

  it('maxMismatches=2 returns at least as many hits as maxMismatches=1', async () => {
    const { sab: sab1, view: view1 } = makeSab(chromNames);
    await searcher.streamAlignments(Buffer.from(sab1), [UNIQUE_QUERY], 1);

    const { sab: sab2, view: view2 } = makeSab(chromNames);
    await searcher.streamAlignments(Buffer.from(sab2), [UNIQUE_QUERY], 2);

    expect(view2.committedCount).toBeGreaterThanOrEqual(view1.committedCount);

    const cursor = view2.allocateCursor();
    for (let seq = 0; seq < view2.committedCount; seq++) {
      cursor.seek(seq);
      expect(cursor.getU8('mismatches')).toBeLessThanOrEqual(2);
    }
  });
});
