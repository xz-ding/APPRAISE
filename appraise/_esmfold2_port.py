"""Start the pinned ColabFold2 port with a bounded-memory ESM-C weight loader."""

import hashlib
import runpy
import sys
from importlib import metadata
from importlib.abc import MetaPathFinder
from importlib.machinery import SourceFileLoader
from importlib.util import spec_from_file_location
from textwrap import indent


ESM_MODULE = "alphafold3.model.esm"
UPSTREAM_SHA256 = "29bbb6954c9840df2d6e7998eaf6ecf08228a804da38329412f85ffc631f2f47"
STREAMING_SHA256 = "304ed05b7717130cdfff74f9f1ce80efcde1c08d3b37b1dc42dddc8a66379efc"


def _streaming_source(data):
    if hashlib.sha256(data).hexdigest() != UPSTREAM_SHA256:
        raise ValueError("The installed ESM-C loader differs from the verified ColabFold2 port.")
    source = data.decode("utf-8")

    # Read the same records directly from Zstandard. Stream one record and free
    # the source tensors for each group as it is stacked.
    start = source.index("  with open(path, 'rb') as fh:\n", source.index("def load("))
    loop = source.index("  for scope, name, arr in read_records(io.BytesIO(raw)):\n", start)
    stack = source.index("  for key, layers in per_block.items():\n", loop)
    streamed = (
        "  with open(path, 'rb') as fh:\n"
        "    stream = zstandard.ZstdDecompressor().stream_reader(fh) if path.endswith('.zst') else fh\n"
        + indent(source[loop:stack].replace("read_records(io.BytesIO(raw))", "read_records(stream)"), "  ")
    )
    source = source[:start] + streamed + source[stack:]
    stack_code = "    p[key] = np.stack([layers[i] for i in sorted(layers)])\n"
    source = source.replace(
        "  for key, layers in per_block.items():\n" + stack_code,
        "  for key in list(per_block):\n    layers = per_block.pop(key)\n" + stack_code + "    del layers\n",
        1,
    )

    # The port writes native NumPy caches. Reopen the freshly written cache as
    # memory maps just as the upstream loader already does on subsequent runs.
    write = "    _write_cache(cache, p)\n"
    source = source.replace(write, write + "    mapped = _read_cache(cache)\n    if mapped is not None:\n      p = mapped\n", 1)
    if hashlib.sha256(source.encode("utf-8")).hexdigest() != STREAMING_SHA256:
        raise ValueError("The streamed ESM-C loader differs from the tested source.")
    return source


class _StreamingLoader(SourceFileLoader):
    def get_code(self, fullname):
        source = _streaming_source(self.get_data(self.path))
        print(f"Using streamed ESM-C weight loader {STREAMING_SHA256}", flush=True)
        return compile(source, self.path, "exec", dont_inherit=True)


class _StreamingPort(MetaPathFinder):
    def find_spec(self, fullname, path=None, target=None):
        if fullname != ESM_MODULE:
            return None
        package = metadata.distribution("alphafold3-colabfold")
        origin = str(package.locate_file("alphafold3/model/esm.py"))
        return spec_from_file_location(fullname, origin, loader=_StreamingLoader(fullname, origin))


if __name__ == "__main__":
    # Apply the loader in this model subprocess. The installed port remains
    # available to other clients through its ordinary import path.
    sys.meta_path.insert(0, _StreamingPort())
    sys.argv = sys.argv[1:]
    runpy.run_path(sys.argv[0], run_name="__main__")
