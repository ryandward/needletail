/**
 * bench.mjs — Performance benchmark: FM-Index build + alignment search.
 *
 * Usage: node bench.mjs <fasta> <query> [query2 ...]
 */

import { buildSchema, computeStrandMap, initStrandHeader, StrandView } from '@strand/core';
import { createRequire } from 'module';
import { statSync } from 'fs';

const require = createRequire(import.meta.url);
const { NativeFmIndex } = require('./index.js');

const schema = buildSchema([
  { name: 'chrom',      type: 'utf8_ref' },
  { name: 'pos',        type: 'u32'      },
  { name: 'query_id',   type: 'u32'      },
  { name: 'strand',     type: 'bool8'    },
  { name: 'mismatches', type: 'u8'       },
  { name: 'score',      type: 'f32'      },
]);

const [,, fastaPath, ...queries] = process.argv;
if (!fastaPath || queries.length === 0) {
  console.error('Usage: node bench.mjs <fasta> <query> [query2 ...]');
  process.exit(1);
}

const genomeSizeMB = (statSync(fastaPath).size / 1_048_576).toFixed(1);
console.log(`Genome : ${fastaPath}  (${genomeSizeMB} MB)`);
console.log(`Queries: ${queries.join(', ')}`);
console.log('');

// ── Index build ───────────────────────────────────────────────────────────────
process.stdout.write('Building FM-Index … ');
const t0 = performance.now();
const searcher = new NativeFmIndex(fastaPath);
const buildMs = (performance.now() - t0).toFixed(0);
console.log(`done  ${buildMs} ms`);

// Chromosome names travel in the SAB header (v4 metadata region) so any
// consumer can resolve handles without a separate side-channel call.
const chromNames = searcher.chromNames();
console.log(`Chroms : ${chromNames.length}  (${chromNames[0]} … ${chromNames.at(-1)})`);
console.log('');

function makeSab() {
  const map = computeStrandMap({
    schema,
    index_capacity: 65536,
    heap_capacity: 0,
    query: { assembly: '', chrom: '', start: 0, end: 0 },
    estimated_records: 0,
  });
  const sab = new SharedArrayBuffer(map.total_bytes);
  // Embed intern table in header — consumers call readStrandHeader(sab).meta.intern
  initStrandHeader(sab, map, { intern: chromNames });
  const view = new StrandView(sab);
  view.updateInternTable(chromNames);
  return { sab, view };
}

// ── Search ────────────────────────────────────────────────────────────────────
const strandLabel = (b) => b ? '+' : '−';

for (const query of queries) {
  const { sab, view } = makeSab();

  const t1 = performance.now();
  await searcher.streamAlignments(Buffer.from(sab), [query]);
  const searchMs = (performance.now() - t1).toFixed(2);

  const committed = view.committedCount;
  console.log(`Query  : ${query}  (${query.length} nt)`);
  console.log(`Search : ${searchMs} ms   Hits: ${committed}`);

  if (committed > 0) {
    const cursor = view.allocateCursor();
    for (let seq = 0; seq < committed; seq++) {
      cursor.seek(seq);
      const chrom = cursor.getRef('chrom');   // NC_001133.9 etc.
      const pos   = cursor.getU32('pos');
      const fwd   = cursor.getBool('strand');
      const mm    = cursor.getU8('mismatches');
      const score = cursor.getF32('score');
      console.log(
        `  [${seq}]  chrom=${chrom}  pos=${pos}  strand=${strandLabel(fwd)}  mm=${mm}  score=${score.toFixed(4)}`
      );
    }
  }
  console.log('');
}
