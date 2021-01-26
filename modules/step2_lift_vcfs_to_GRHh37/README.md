# Lift vcfs to GRHh37

This module takes a vcf file as input and uses the build information to lift the file to Build37 using Plink .


## Testing

Since DSL2 allows a module to have a workflow as well, we can test this module and optimize directives using the `test` workflow.

This process could be tested using the `test` workflow which relies on

- the local files in `test_data`
- params in  `test_params.yaml`
- `test` profile in `nextflow.config`


```
nextflow run lift_vcfs_to_GCRHh37.nf -entry test -params-file test_params.yaml -entry test
```

We can use the optimum configuration without making any change to the module file, for the production workloads, by simply copy and paste the desired configs to the workflow file or workflow params.


- Workflow file e.g. `nextflow.config`

<!-- ```nextflow
params.CHECK_VCF_BUILD [
// add the optimized params
java_opts : ""
]
include { CHECK_VCF_BUILD } from "relative_path_to/modules/step1-check_vcf_build/check_vcf_build.nf" addParams (params.CHECK_VCF_BUILD)
``` -->

- Workflow params file e.g. `params.yaml`

<!-- ```yaml

CHECK_VCF_BUILD:
  java_opts : ""

``` -->
