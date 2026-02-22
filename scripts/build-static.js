const fs = require("node:fs/promises");
const path = require("node:path");

const ROOT = process.cwd();
const PUBLIC_DIR = path.join(ROOT, "public");
const DIST_DIR = path.join(ROOT, "dist");

async function main() {
  // Recreate dist/
  await fs.rm(DIST_DIR, { recursive: true, force: true });
  await fs.mkdir(DIST_DIR, { recursive: true });

  // Copy public/* -> dist/*
  await fs.cp(PUBLIC_DIR, DIST_DIR, { recursive: true });

  // Basic sanity check
  await fs.access(path.join(DIST_DIR, "index.html"));
  console.log(`Built static site: ${path.relative(ROOT, DIST_DIR)}/`);
}

main().catch((err) => {
  console.error(err);
  process.exit(1);
});


