/**
 * align.mjs — Live alignment demonstration against phiX174 genome.
 *
 * Usage: node align.mjs [query1] [query2] ...
 * Default queries are taken from known phiX174 subsequences.
 */

import { buildSchema, computeStrandMap, initStrandHeader, StrandView } from '@strand/core';
import { createRequire } from 'module';
import { resolve, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const require = createRequire(import.meta.url);

const { NativeFmIndex } = require('./index.js');

// ── Schema — must match strand.rs constants exactly ───────────────────────────
const schema = buildSchema([
  { name: 'chrom_id',   type: 'u32'   },
  { name: 'pos',        type: 'u32'   },
  { name: 'strand',     type: 'bool8' },
  { name: 'mismatches', type: 'u8'    },
  { name: 'score',      type: 'f32'   },
]);

function makeSab(indexCapacity = 65536) {
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

// ── Queries from phiX174 (NC_001422.1, 5386 nt) ───────────────────────────────
// These are real substrings extracted from the phiX174 genome.
const DEFAULT_QUERIES = [
  // Forward strand, pos 0: canonical phiX174 start of gene A
  'GAGTTTTATCGCTTCCATGACGCAG',
  // Forward strand, pos 115: within gene A
  'TGGACTGCTGGCGGAAAATGAGAAA',
  // Forward strand, pos 389: within gene A
  'ATGAGTCAAGTT',
  // Minus strand query: RC of seq[2000:2020] → should align at pos=2000, strand=-
  'CGAATCACCAGAACGGAAAA',
  // Absent control
  'NNNNNNNNNNNNNNNNNNNN',
];

const queries = process.argv.length > 2
  ? process.argv.slice(2)
  : DEFAULT_QUERIES;

const GENOME = resolve(__dirname, 'test/fixtures/phiX174.fa');

// ── Run ───────────────────────────────────────────────────────────────────────
console.log(`Genome : ${GENOME}`);
console.log(`Queries: ${queries.length}`);
console.log('');

const searcher = new NativeFmIndex(GENOME);

let totalHits = 0;

for (const query of queries) {
  const { sab } = makeSab();
  const view = new StrandView(sab);

  const t0 = performance.now();
  await searcher.streamAlignments(Buffer.from(sab), [query]);
  const elapsed = (performance.now() - t0).toFixed(2);

  const committed = view.committedCount;
  totalHits += committed;

  const strandLabel = (b) => b ? '+' : '−';

  if (committed === 0) {
    console.log(`Query : ${query}`);
    console.log(`Result: 0 hits  (${elapsed} ms)`);
  } else {
    console.log(`Query : ${query}`);
    console.log(`Hits  : ${committed}  (${elapsed} ms)`);

    const cursor = view.allocateCursor();
    for (let seq = 0; seq < committed; seq++) {
      cursor.seek(seq);
      const chrom = cursor.getU32('chrom_id');
      const pos   = cursor.getU32('pos');
      const fwd   = cursor.getBool('strand');
      const mm    = cursor.getU8('mismatches');
      const score = cursor.getF32('score');
      console.log(
        `  [${seq}] chrom=${chrom}  pos=${pos}  strand=${strandLabel(fwd)}  mm=${mm}  score=${score.toFixed(4)}`
      );
    }
  }
  console.log('');
}

console.log(`Total hits across all queries: ${totalHits}`);
