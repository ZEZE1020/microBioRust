 <img src="docs/assets/MICROBIO B.svg" width=200em alt="logo" /> 

[![Docs](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://lcrossman.github.io/microBioRust/)

![Crates.io Version](https://img.shields.io/crates/v/microBioRust?style=flat&link=https%3A%2F%2Fcrates.io%2Fcrates%2FmicroBioRust)




## A Rust bioinformatics crate aimed at Microbial genomics<br>

The aim of this crate is to provide Microbiology friendly Rust functions for bioinformatics.<br>

> Very much under construction!<br>

Some concepts with many thanks to Rust-bio<br>
Please see the Roadmap for futher details [here](ROADMAP.md)

To install Rust - please see here [Rust install](https://www.rust-lang.org/tools/install) or with Conda<br>
If you would like to contribute please follow the [Rust code of conduct](https://www.rust-lang.org/policies/code-of-conduct)

Questions and comments - please join the Discord server :) [here](https://discord.gg/xP2ngwTttz)


Currently there is functionality for:<br>
````
 1. A Genbank to GFF parser

 2. An Embl to GFF and GBK parser

 3. Calculate sequence metrics e.g. hydrophobicity, distance measures

 4. A Heatmap plot with wasm and d3.js

````

To see more on how to use have a look at usage [here](docs/usage.md)

To use a specific workspace (at the moment microSeqIO or heatmap) clone the project, cd into the specific directory required and build the project from there

For more background please see <https://LCrossman.github.io/microBioRust_details>
