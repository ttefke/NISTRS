use libm::lgamma;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use statrs::function::gamma::gamma_ur;

mod template10;
mod template11;
mod template12;
mod template13;
mod template14;
mod template15;
mod template16;
mod template2;
mod template3;
mod template4;
mod template5;
mod template6;
mod template7;
mod template8;
mod template9;

use self::{
    template10::TEMPLATE10, template11::TEMPLATE11, template12::TEMPLATE12, template13::TEMPLATE13,
    template14::TEMPLATE14, template15::TEMPLATE15, template16::TEMPLATE16, template2::TEMPLATE2,
    template3::TEMPLATE3, template4::TEMPLATE4, template5::TEMPLATE5, template6::TEMPLATE6,
    template7::TEMPLATE7, template8::TEMPLATE8, template9::TEMPLATE9,
};

use super::*;

const MAX_NUM_OF_TEMPLATES: usize = 47914;

/// Overlapping Template Matching Test.
/// The focus of the Overlapping Template Matching test is the number of occurrences of pre-specified target
/// strings. Both this test and the Non-overlapping Template Matching test use an m-bit
/// window to search for a specific m-bit pattern. If the pattern is not found,
/// the window slides one bit position. The difference between this test and the test in Non-overlapping Template Matching test is that
/// when the pattern is found, the window slides only one bit before resuming the search.
/// `m` - the length in bits of each template.
pub fn overlapping_template_test(data: &BitsData, m: usize) -> Result<Vec<TestResultT>, String> {
    const M: usize = 1032;
    const N: usize = 968;

    let n = data.len() / M;

    if data.len() < M * N {
        return Err(format!("n is too small, increase to at least {}", M * 968).to_string());
    }

    let pi = compute_pi(M, m);

    let num_of_rows = get_template_size(m).min(MAX_NUM_OF_TEMPLATES);
    let par_iter = (0..num_of_rows).into_par_iter().map(|i| {
        let test_seq = get_template(i, m);

        let mut nu: [f64; 6] = [0_f64; 6];
        let mut begin = usize::default();
        let mut end = begin + m;

        for _ in 0..n {
            let mut w_obs = usize::default();
            for _ in 0..(M - m + 1) {
                if data[begin..end].eq(test_seq) {
                    w_obs += 1;
                }
                begin += 1;
                end += 1;
            }

            begin += m - 1;
            end += m - 1;
            w_obs = w_obs.min(nu.len() - 1);
            nu[w_obs] += 1_f64;
        }
        let mut chi2 = 0_f64;
        nu.iter().zip(pi.iter()).for_each(|x| {
            chi2 += (x.0 - x.1 * (n as f64)).powi(2) / (x.1 * (n as f64));
        });

        let p = gamma_ur(2.5_f64, chi2 / 2_f64);

        (p >= TEST_THRESHOLD, p)
    });

    Ok(par_iter.collect())
}

fn compute_pi(c_m: usize, m: usize) -> [f64; 6] {
    let lambda = (c_m - m + 1) as f64 / 2_f64.powi(m as i32);
    let eta = lambda / 2_f64;

    let mut sum = 0_f64;
    let mut res: [f64; 6] = [0_f64; 6];

    for (i, it) in res.iter_mut().enumerate().take(5) {
        let tmp = single_compute_pi(i, eta);
        sum += tmp;

        *it = tmp;
    }

    res[5] = 1_f64 - sum;

    res
}

fn single_compute_pi(u: usize, eta: f64) -> f64 {
    if u == 0 {
        return (-eta).exp();
    }

    let mut sum = 0_f64;
    for l in 1..=u {
        sum += (-eta - (u as f64) * 2_f64.ln() + (l as f64) * eta.ln() - lgamma(l as f64 + 1_f64)
            + lgamma(u as f64)
            - lgamma(l as f64)
            - lgamma((u as f64) - (l as f64) + 1_f64))
        .exp();
    }

    sum
}

fn get_template(i: usize, m: usize) -> &'static [bool] {
    match m {
        2 => &TEMPLATE2[i],
        3 => &TEMPLATE3[i],
        4 => &TEMPLATE4[i],
        5 => &TEMPLATE5[i],
        6 => &TEMPLATE6[i],
        7 => &TEMPLATE7[i],
        8 => &TEMPLATE8[i],
        9 => &TEMPLATE9[i],
        10 => &TEMPLATE10[i],
        11 => &TEMPLATE11[i],
        12 => &TEMPLATE12[i],
        13 => &TEMPLATE13[i],
        14 => &TEMPLATE14[i],
        15 => &TEMPLATE15[i],
        16 => &TEMPLATE16[i],
        _ => panic!("Unknown size! M = {}", m),
    }
}

fn get_template_size(m: usize) -> usize {
    match m {
        2 => TEMPLATE2.len(),
        3 => TEMPLATE3.len(),
        4 => TEMPLATE4.len(),
        5 => TEMPLATE5.len(),
        6 => TEMPLATE6.len(),
        7 => TEMPLATE7.len(),
        8 => TEMPLATE8.len(),
        9 => TEMPLATE9.len(),
        10 => TEMPLATE10.len(),
        11 => TEMPLATE11.len(),
        12 => TEMPLATE12.len(),
        13 => TEMPLATE13.len(),
        14 => TEMPLATE14.len(),
        15 => TEMPLATE15.len(),
        16 => TEMPLATE16.len(),
        _ => panic!("Unknown size! M = {}", m),
    }
}
