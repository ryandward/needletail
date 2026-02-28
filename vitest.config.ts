import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    environment: 'node',
    // Allow SharedArrayBuffer (requires cross-origin isolation in browsers,
    // but Node.js enables it unconditionally).
    pool: 'forks',
    poolOptions: {
      forks: {
        // Each test file gets its own process to avoid SAB cross-contamination.
        singleFork: false,
      },
    },
  },
});
