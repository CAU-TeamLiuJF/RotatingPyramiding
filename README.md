# RotatingPyramiding

----

### Simulate a strategy

```shell
Rscript run.R -p $param_file -r $repeat -s $strategy_name -o $output
```

`$param_file` is the param file name in directory `./params`, default value is `example.txt`.

`$repeat` is the repeat times of simulation, default value is `20`.

`$strategy_name` is the strategy R script name in directory `./strategy` without `.R`, default value is `baseline`.

`$output` is the tmp data directory in directory `./tmp`, default value is `default_output`.

### Get results of a strategy

```shell
Rscript get_result.R -r $repeat -s $strategy_name -o $output
```

`$output` is the output result directory in directory `./output`, default value is `default_output`.

The directory `./tmp` and `./output` will be created automatically. So does `$output`.

### Get target gene distribution results of a strategy

```shell
Rscript get_target_gene.R -r $repeat -s $strategy_name -o $output
```

