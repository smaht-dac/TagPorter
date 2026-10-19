use criterion::{criterion_group, criterion_main, Criterion};

fn my_function() {
    // Code to benchmark
    let _: Vec<_> = (0..1000).map(|x| x * 2).collect();
}

fn bench_my_function(c: &mut Criterion) {
    c.bench_function("my_function", |b| b.iter(|| my_function()));
}

fn prime_factors(mut n: u64) -> Vec<u64> {
    let mut factors = Vec::new();
    let mut d = 2;

    while d * d <= n {
        while n % d == 0 {
            factors.push(d);
            n /= d;
        }
        d += 1;
    }

    if n > 1 {
        factors.push(n);
    }

    factors
}


fn bench_multiple_inputs(c: &mut Criterion) {
    let mut group = c.benchmark_group("vec sizes");

    for size in [100, 200, 300, 400] {
        group.bench_with_input(format!("size={}", size), &size, |b, &s| {
            b.iter(|| {
                let _ = prime_factors(s);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_my_function, bench_multiple_inputs);
criterion_main!(benches);
