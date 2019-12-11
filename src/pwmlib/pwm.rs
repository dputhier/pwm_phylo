extern crate regex;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufRead;
use std::io::BufReader;

use regex::Regex;

pub fn read_pwm_as_tab(input_file: &str,
                       pwm_to_score: &mut HashMap<(char, i32), f64>) -> i32 {


    // pwm size
    let mut pwm_size = 0;

    // open the file
    self::super::utils::message(&format!("Found position-weight matrix: {}.", input_file));

    let _file_handler = File::open(input_file).unwrap();


    // Loop over lines
    let _buffer_reader = BufReader::new(&_file_handler);

    for (_num, _result) in _buffer_reader.lines().enumerate() {
        let line = _result.unwrap();

        // get the nucleotide name
        let _nuc = line.chars().next().unwrap();

        // replace the beginning of the string
        // Recall the tab format:
        //a       |  5 4 15...
        //c       |  3 3 3...
        //...

        let re_rep = Regex::new(r"^.*\|\s*").unwrap();
        let line = re_rep.replace_all(&line, "");

        // Split the line using ' ' or '  '
        let re_split = Regex::new(r"\s+").unwrap();

        // Loop over the splitted line
        pwm_size = 0;

        for val_str in re_split.split(&line) {
            let _val_num: f64 = val_str.parse().unwrap();
            pwm_to_score.insert((_nuc, pwm_size), _val_num);
            pwm_size += 1;
        }
    }


    pwm_size
}
