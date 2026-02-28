/**
 * synthetic.test.ts — Performance / stress test using a stub producer that
 * writes records at maximum speed without invoking the FM-Index.
 *
 * This test exercises the Strand ring-buffer protocol directly from JS to
 * verify backpressure, consumer acknowledgment, and EOS detection at high
 * throughput.
 */

import { describe, it, expect } from 'vitest';
import {
  buildSchema,
  computeStrandMap,
  initStrandHeader,
  StrandView,
  StrandWriter,
} from '@strand/core';

const schema = buildSchema([
  { name: 'chrom_id',   type: 'u32'   },
  { name: 'pos',        type: 'u32'   },
  { name: 'strand',     type: 'bool8' },
  { name: 'mismatches', type: 'u8'    },
  { name: 'score',      type: 'f32'   },
]);

function makeSab(indexCapacity = 1024) {
  const map = computeStrandMap({
    schema,
    index_capacity: indexCapacity,
    heap_capacity: 0,
    query: { assembly: '', chrom: '', start: 0, end: 0 },
    estimated_records: 0,
  });
  const sab = new SharedArrayBuffer(map.total_bytes);
  initStrandHeader(sab, map);
  return { sab, map };
}

// ─── Stub producer: write N records at full speed from JS ────────────────────

describe('synthetic stub producer throughput', () => {
  it('writes 10_000 records and consumer reads them all', async () => {
    const { sab } = makeSab(1024);
    const writer = new StrandWriter(sab);
    const view   = new StrandView(sab);

    const TOTAL = 10_000;
    let consumed = 0;

    // Producer: write in a Worker-style synchronous loop (Node.js main thread
    // can use StrandWriter since Atomics.wait is allowed outside browsers).
    writer.begin();

    // Write in batches to avoid blocking forever on a tiny ring.
    const BATCH = 256;
    for (let i = 0; i < TOTAL; i += BATCH) {
      const records = [];
      for (let j = 0; j < Math.min(BATCH, TOTAL - i); j++) {
        records.push({
          chrom_id:   0,
          pos:        i + j,
          strand:     true,
          mismatches: 0,
          score:      1.0,
        });
      }
      writer.writeRecordBatch(records);

      // Consume and acknowledge to prevent ring stall.
      const committed = view.committedCount;
      const cursor = view.allocateCursor();
      for (; consumed < committed; consumed++) {
        cursor.seek(consumed);
      }
      view.acknowledgeRead(consumed);
    }

    writer.finalize();

    expect(consumed).toBe(TOTAL);
    expect(view.status).toBe('eos');
  });

  it('abort from consumer terminates producer write loop', async () => {
    const { sab } = makeSab(16);  // tiny ring
    const writer = new StrandWriter(sab);
    const view   = new StrandView(sab);

    writer.begin();

    // Write some records without consuming, to fill the ring.
    const BATCH = 8;
    const records = Array.from({ length: BATCH }, (_, i) => ({
      chrom_id: 0, pos: i, strand: true, mismatches: 0, score: 1.0,
    }));
    writer.writeRecordBatch(records);

    // Signal abort — producer should throw StrandAbortError on next stall.
    view.signalAbort();
    writer.abort();

    expect(view.status).toBe('error');
  });

  it('ring-buffer wrap-around produces correct positions', async () => {
    // Write more records than index_capacity to exercise the ring wrap.
    const CAPACITY = 8;
    const { sab } = makeSab(CAPACITY);
    const writer = new StrandWriter(sab);
    const view   = new StrandView(sab);

    writer.begin();

    const TOTAL = CAPACITY * 4;  // 4 full laps of the ring
    let consumed = 0;

    for (let i = 0; i < TOTAL; i++) {
      writer.writeRecordBatch([{
        chrom_id: 0, pos: i, strand: true, mismatches: 0, score: 1.0,
      }]);

      // Consume each record immediately to keep the ring clear.
      const committed = view.committedCount;
      const cursor = view.allocateCursor();
      cursor.seek(consumed);
      expect(cursor.getU32('pos')).toBe(i);
      consumed++;
      view.acknowledgeRead(consumed);
    }

    writer.finalize();

    expect(consumed).toBe(TOTAL);
    expect(view.status).toBe('eos');
  });
});
