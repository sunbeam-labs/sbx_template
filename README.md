# Sunbeam extension template

This is a template to use to extend the [Sunbeam pipeline](https://github.com/sunbeam-labs/sunbeam). There are three major parts to a Sunbeam extension: 

 - `requirements.txt` specifies the extension's dependencies
 - `config.yml` contains configuration options that can be specified by the user when running an extension
 - `sbx_template.rules` contains the rules (logic/commands run) of the extension
 
## Creating an extension

Any dependencies required for an extension should be listed in the `requirements.txt` file in the [standard requirement format](https://pip.readthedocs.io/en/1.1/requirements.html). 

The `config.yml` contains parameters that the user might need to modify when running an extension. For example, if your downstream analysis is run differently depending on whether reads are paired- or single-end, it would probably be wise to include a `paired_end` parameter. Default values should be specified for each bottom-level key.

Finally, `sbx_template.rules` contains the actual logic for the extension, including required input and output files. A detailed discussion of Snakemake rule creation is beyond the scope of this tutorial, but definitely check out [the Snakemake tutorial](http://snakemake.readthedocs.io/en/stable/tutorial/basics.html) and any of the [extensions by sunbeam-labs](https://github.com/sunbeam-labs) for inspiration.

## Installing an extension

Installing an extension is as simple as cloning (or moving) your extension directory into the sunbeam/extensions/ folder, installing requirements through Conda, and adding the new options to your existing configuration file: 

    git clone https://github.com/louiejtaylor/sbx_template/ sunbeam/extensions/sbx_template
    conda install --file sunbeam/extensions/sbx_template/requirements.txt
    cat sunbeam/extensions/sbx_template/config.yml >> sunbeam_config.yml

## Running an extension

To run an extension, simply run Sunbeam as usual with your extension's target rule specified:

    sunbeam run --configfile=sunbeam_config.yml example_rule
    
