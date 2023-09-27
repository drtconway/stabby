# Stabby

This library implements Schidt's efficient data structure for integer stabbing in Rust.

https://doi.org/10.1007/978-3-642-10631-6_18

The basic usage is as follows:

```rust
use stabby::{Interval, Stabby};

// Make some random intervals:
//
use rand::Rng;
let mut rng = rand::thread_rng();
let mut intervals = Vec::new();
for i in 0..100 {
    let first = rng.gen_range(1, 10000);
    let last = first + rng.gen_range(1, 500);
    intervals.push(Interval::new(first, last));
}
intervals.sort();

// Now build the stabby index:
//
let idx = Stabby::new(&intervals);

// And stab a single position
//
let query = rng.gen_range(1, 10000);
let results = idx.stab(query);
for result in results.iter() {
    //check for overlap.
    assert!(result.first <= query && query <= result.last);
}

// Or stab an interval
//
let query_first = rng.gen_range(1, 10000);
let query_last = query_first + rng.gen_range(1, 500);
let query = Interval::new(query_first, query_last);
let results = idx.stab_interval(&query);
for result in results.iter() {
    // check for overlap.
    assert!(query.first <= result.last && result.first <= query.last);
}
```
