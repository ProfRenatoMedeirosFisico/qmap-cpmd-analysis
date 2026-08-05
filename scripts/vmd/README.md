# VMD/Tcl scripts

Place here the exact `.tcl` scripts used in the published analysis.

Each script should include:

- required VMD version;
- trajectory/topology loading command;
- atom selections;
- periodic-boundary handling;
- frame range and stride;
- units;
- output filename;
- command used to execute the script, for example:

```bash
vmd -dispdev text -e scripts/vmd/analysis_name.tcl
```

Do not replace the original analysis scripts with newly reconstructed scripts without clearly documenting and validating the change.
