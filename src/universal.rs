use libm::erfc;

use super::*;

/// Maurer’s “Universal Statistical”.
/// The focus of this test is the number of bits between matching patterns (a measure that is related to the
/// length of a compressed sequence). The purpose of the test is to detect whether or not the sequence can be
/// significantly compressed without loss of information. A significantly compressible sequence is
/// considered to be non-random.
pub fn universal_test(data: &BitsData) -> Result<[TestResultT; 4], String> {
    const EXPECTED_VALUE: [f64; 17] = [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 5.2177052, 6.1962507, 7.1836656, 8.1764248, 9.1723243,
        10.170032, 11.168765, 12.168070, 13.167693, 14.167488, 15.167379,
    ];

    const VARIANCE: [f64; 17] = [
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384, 3.401, 3.410,
        3.416, 3.419, 3.421,
    ];

    let n_bits = data.len();

    let l: usize = match n_bits {
        1059061760.. => 16,
        496435200.. => 15,
        231669760.. => 14,
        107560960.. => 13,
        49643520.. => 12,
        22753280.. => 11,
        10342400.. => 10,
        4654080.. => 9,
        2068480.. => 8,
        904960.. => 7,
        387840.. => 6,
        _ => 5,
    };

    // const double c = 0.7 - 0.8/static_cast<double>(L)
    //         + (4.0 + 32.0/static_cast<double>(L)) * std::pow(K, -3.0/static_cast<double>(L))/15.0;

    let p = 2_usize.pow(l as u32);
    let q = 10 * p;
    let k = n_bits / l - q;

    if n_bits < (q + k) * l {
        return Err("n is too small!".to_string());
    }

    if l < 6 || l > 16 {
        return Err("Invalid value of l, check size of n!".to_string());
    }

    let c = [
        0.7 - 0.8 / (l as f64)
        + (4_f64 + 32_f64 / (l as f64)) * (k as f64).powf(-3_f64 / (l as f64)) / 15_f64,
        0.7 - 0.8 / (l as f64)
        + (1.6_f64 + 12.8_f64 / (l as f64)) * (k as f64).powf(-4_f64 / (l as f64))
    ];

    let sigma = [
        c[0] * (VARIANCE[l] / (k as f64)).sqrt(),
        c[1] * (VARIANCE[l] / (k as f64)).sqrt()
    ];

    let mut t = vec![usize::default(); p];

    for i in 1..=q {
        let mut dec_rep = usize::default();
        let ind_i = i - 1;
        for j in 0..l {
            if data[ind_i * l + j] {
                dec_rep += 2_usize.pow((l - 1 - j) as u32);
            }
        }

        t[dec_rep] = i;
    }

    let mut sum = f64::default();
    for i in q + 1..=q + k {
        let mut dec_rep = usize::default();
        let ind_i = i - 1;
        for j in 0..l {
            if data[ind_i * l + j] {
                dec_rep += 2_usize.pow((l - 1 - j) as u32);
            }
        }
        sum += ((i - t[dec_rep]) as f64).ln() / 2_f64.ln();
        t[dec_rep] = i;
    }

    let phi = sum / (k as f64);
    let arg = [
        (phi - EXPECTED_VALUE[l]).abs() / (2_f64.sqrt() * sigma[0]),
        (phi - EXPECTED_VALUE[l]).abs() / (2_f64.sqrt() * sigma[1]),
        // t-test
        (phi - EXPECTED_VALUE[l]).abs() / sigma[0],
        (phi - EXPECTED_VALUE[l]).abs() / sigma[1],
    ];

    let p = [
        erfc(arg[0]),
        erfc(arg[1]),
        erfc(arg[2]),
        erfc(arg[3])
    ];

    Ok([
        (p[0] >= TEST_THRESHOLD, p[0]),
        (p[1] >= TEST_THRESHOLD, p[1]),
        (p[2] >= TEST_THRESHOLD, p[2]),
        (p[3] >= TEST_THRESHOLD, p[3])
    ])
}
