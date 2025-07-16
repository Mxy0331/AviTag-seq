# AviTag-seq

AviTag-seq is a collection of scripts used to process AviTag sequencing data and
identify potential off-target integration sites.

## Requirements

Python 3.8 or above is recommended. External tools such as **bwa** and
**bedtools** are used during alignment steps. Python dependencies are listed in
[`requirements.txt`](requirements.txt).

Install the Python dependencies using pip:

```bash
pip install -r requirements.txt
```

## Usage

The main workflow is implemented in `AviTag-Seq_0411/AviTag-seq.py` and provides
a command line interface. The typical command to run all steps is:

```bash
python AviTag-Seq_0411/AviTag-seq.py all --manifest path/to/manifest.yaml
```

To run only the alignment and analysis steps without read splitting you can
use the `findoff` subcommand:

```bash
python AviTag-Seq_0411/AviTag-seq.py findoff --manifest path/to/manifest.yaml
```

See the manifest example in `AviTag-Seq_0411/manifest.yaml` for the required
fields.

## Tests

Basic unit tests can be executed with `pytest`:

```bash
pytest
```


