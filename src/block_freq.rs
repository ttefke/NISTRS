use statrs::function::gamma::gamma_ur;

use super::*;

/// Frequency Test within a Block.
///
/// The focus of the test is the proportion of ones within M-bit blocks. The purpose of this test is to determine
/// whether the frequency of ones in an M-bit block is approximately M/2, as would be expected under an
/// assumption of randomness. For block size M=1, this test degenerates to test 1, the Frequency (Monobit) test.
/// `M` must be > 0;
/// # Example
/// ```
/// use nistrs::block_freq::block_frequency_test;
/// use nistrs::BitsData;

/// let data = BitsData::from_text("11001001000011111101101010100010001000010110100011000010001101001100010011000110011000101000101110000000".to_string());
/// assert_eq!(block_frequency_test(&data, 10).unwrap().1, 0.70643844964128211);
/// ```
pub fn block_frequency_test(data: &super::BitsData, m: usize) -> Result<TestResultT, String> {
    let nbits = data.len();

    if nbits < 100 {
        return Err("At least 100 bits are required to run the block frequency test.".to_string());
    } else if nbits < m {
        return Err("Size of block must be lower than number of bits".to_string());
    }

    let n_blocks = nbits / m;

    if nbits < n_blocks * m {
        return Err(
            "Number of provided bits is smaller than amount of storage in the blocks.".to_string(),
        );
    }

    if m < 20 {
        return Err("The block size m should be at least 20.".to_string());
    }

    let min_m = 0.01 * nbits as f64;
    if m as f64 <= min_m {
        return Err(format!(
            "Each block must hold more than 1% of the input, increase m to at least {}",
            min_m as u64 + 1
        ));
    }

    if n_blocks >= 100 {
        return Err(format!(
            "Too many blocks, increase m to a maximum of {}",
            nbits / 99
        ));
    }

    let mut sum = f64::default();
    for i in 0..n_blocks {
        let mut block_sum = 0_usize;
        for j in 0..m {
            block_sum += data[i * m + j] as usize;
        }

        let v = (block_sum as f64) / (m as f64) - 0.5;
        sum += v.powf(2_f64);
    }

    let chi_squared = 4_f64 * (m as f64) * sum;
    let p = gamma_ur((n_blocks as f64) / 2_f64, chi_squared / 2_f64);

    Ok((p >= TEST_THRESHOLD, p))
}
