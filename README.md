# Great Circle Triple Finder

Simple Node app with two views:
1. A cutoff-filtered, scrollable table of deduplicated location sets:
   - Build great circles from every pair of locations.
   - For each pair-circle, include all locations within the mile cutoff (default 10 miles).
   - Keep only sets with at least `m` locations, where `m` is configurable (min/default 3).
   - Deduplicate sets by location membership (order-insensitive).
   - Show each set with how many pairs generate it and one example pair.
   - Plot each set's best-fit great circle on the map.
   - There is also a checkbox to plot all individual triplet circles currently below the cutoff.
2. A multi-triplet plotter where you can pick several triplets, draw their great circles, and
   list all pairwise circle intersections as latitude/longitude coordinates.

## Input format

One location per line:

`<place name>, LAT, LONG`

`LAT/LONG` can be plain signed decimals (e.g. `36.7538, 3.0588`) or hemisphere style
(e.g. `36.7538° N, 3.0588° E`).

Example:

`San Francisco, 37.7749, -122.4194`

## Run

```bash
npm start
```

Then open `http://localhost:3000`.

## Alternating word-circle algorithm

Use the **Alternating Word Great Circles** panel in the app UI, or call
`buildAlternatingGreatCircles(set1, set2, options)` in the browser console.

- `set1` / `set2`: arrays of words or newline/semicolon-separated strings
- `options.pairsPerCircle` (optional): minimum number of `(set1,set2)` pairs per circle (default `3`)
- Rules enforced:
  - all Set 2 words are used exactly once
  - each word is used at most once
  - words alternate in each circle as `set1, set2, set1, set2, ...`
  - if entries include `<name>, LAT, LONG`, circles are optimized to reduce distance to the
    best-fit great circle

Example:

```js
buildAlternatingGreatCircles(
  ["alpha", "beta", "gamma", "delta", "epsilon"],
  ["one", "two", "three", "four"],
  { pairsPerCircle: 2 }
);
```
