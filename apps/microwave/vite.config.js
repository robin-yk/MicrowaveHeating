import { defineConfig } from "vite";
import { readFileSync } from "node:fs";

// Keep exported assets relative so the app can be hosted under a repository path.
export default defineConfig({
  base: "./",
  plugins: [{
    name: "include-citation-download",
    generateBundle() {
      this.emitFile({
        type: "asset",
        fileName: "citation.ris",
        source: readFileSync(new URL("./citation.ris", import.meta.url)),
      });
    },
  }],
});
