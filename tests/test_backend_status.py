"""Why a compiled backend can or cannot be used: the build record
(buildinfo.py), the states of cppext.py and the report of doctor.py."""
import os
import sysconfig
import warnings

import pytest

from dmrgpy import buildinfo, cppext, doctor


def _backend_folder(tmp_path):
    (tmp_path / "bindings.cc").write_text("// bindings\n")
    (tmp_path / "chain_session.h").write_text("// session\n")
    (tmp_path / "extra").mkdir()
    (tmp_path / "extra" / "all.h").write_text("// extra\n")
    return str(tmp_path)


def test_sources_hash_sees_every_source_and_nothing_else(tmp_path):
    d = _backend_folder(tmp_path)
    h0 = buildinfo.sources_hash(d)
    (tmp_path / "notes.txt").write_text("not a source")
    assert buildinfo.sources_hash(d) == h0
    (tmp_path / "extra" / "all.h").write_text("// extra, edited\n")
    h1 = buildinfo.sources_hash(d)
    assert h1 != h0
    (tmp_path / "new.h").write_text("")
    assert buildinfo.sources_hash(d) != h1


def test_no_sources_means_no_hash(tmp_path):
    assert buildinfo.sources_hash(str(tmp_path)) is None
    assert buildinfo.is_stale(str(tmp_path)) is None


def test_stamp_keeps_the_other_fields_and_stale_follows_the_sources(tmp_path):
    d = _backend_folder(tmp_path)
    assert buildinfo.is_stale(d) is None  # no record yet
    buildinfo.stamp(d, compiler="/usr/bin/g++", modules=["gcc/12"])
    buildinfo.stamp(d)  # what "make pybind" does after a later rebuild
    info = buildinfo.read(d)
    assert info["compiler"] == "/usr/bin/g++" and info["modules"] == ["gcc/12"]
    assert buildinfo.is_stale(d) is False
    (tmp_path / "chain_session.h").write_text("// session, edited\n")
    assert buildinfo.is_stale(d) is True
    buildinfo.stamp(d)
    assert buildinfo.is_stale(d) is False


def test_a_broken_record_reads_as_no_record(tmp_path):
    d = _backend_folder(tmp_path)
    (tmp_path / buildinfo.INFO_NAME).write_text("{ not json")
    assert buildinfo.read(d) is None
    assert buildinfo.is_stale(d) is None


@pytest.fixture
def fake_backend(tmp_path, monkeypatch):
    """Point cppext at an empty folder for the compiled backends and give
    it a clean memory, restored afterwards."""
    monkeypatch.setattr(cppext, "extension_dir", lambda version: str(tmp_path))
    monkeypatch.setattr(cppext, "_states", {})
    monkeypatch.setattr(cppext, "_backends", {3: None, 2: None})
    return tmp_path


def test_no_extension_file_is_not_compiled_and_is_quiet(fake_backend):
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cppext._record_failure(3, ImportError("cannot import name"))
    assert cppext.status(3)["state"] == cppext.NOT_COMPILED
    assert "is not compiled" in cppext.explain(3)
    # a pip install has no install.py: the sentence must say where it is
    assert "clone of the dmrgpy repository" in cppext.explain(3)


def test_an_extension_for_another_python_is_named_and_warned(fake_backend):
    (fake_backend / "_dmrgcpp.cpython-27-nowhere.so").write_text("")
    with pytest.warns(RuntimeWarning, match="compiled for another Python"):
        cppext._record_failure(3, ImportError("cannot import name"))
    assert cppext.status(3)["state"] == cppext.OTHER_PYTHON
    assert "_dmrgcpp.cpython-27-nowhere.so" in cppext.explain(3)


def test_an_extension_that_does_not_load_says_why(fake_backend, monkeypatch):
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    (fake_backend / ("_dmrgcpp"+suffix)).write_text("")
    buildinfo.write(str(fake_backend), dict(
            modules=["gcc/12.3.0", "openblas/0.3"], conda_prefix="/envs/dm"))
    monkeypatch.setenv("LOADEDMODULES", "openblas/0.3:other/1")
    error = ImportError("libstdc++.so.6: version `CXXABI_1.3.15' not found")
    with pytest.warns(RuntimeWarning, match="CXXABI_1.3.15") as record:
        cppext._record_failure(3, error)
    text = str(record[0].message)
    # only what is missing now is named: the module that is not loaded and
    # the conda environment that is not the active one
    assert "gcc/12.3.0" in text and "openblas/0.3" not in text
    assert "/envs/dm" in text
    # the default backend also says what runs meanwhile
    assert "pure-Python backend" in text
    st = cppext.status(3)
    assert st["state"] == cppext.LOAD_FAILED and not st["available"]
    assert st["build"]["modules"] == ["gcc/12.3.0", "openblas/0.3"]


def test_a_failed_first_import_does_not_import_again(fake_backend,
        monkeypatch):
    """The warning is written while the backend is still being loaded.
    Asking for the backend there must not start the import again: on a
    real extension for another Python this recursed without end."""
    monkeypatch.setattr(cppext, "_backends", {})  # nothing loaded yet
    calls = []
    monkeypatch.setattr(cppext, "_load",
            lambda version: calls.append(version))
    (fake_backend / "_dmrgcpp.cpython-27-nowhere.so").write_text("")
    with pytest.warns(RuntimeWarning, match="another Python"):
        cppext._record_failure(3, ImportError("cannot import name"))
    assert calls == [] and cppext._backends[3] is None


def test_v2_warning_does_not_talk_about_the_default(fake_backend):
    suffix = sysconfig.get_config_var("EXT_SUFFIX")
    (fake_backend / ("_dmrgcpp"+suffix)).write_text("")
    with pytest.warns(RuntimeWarning) as record:
        cppext._record_failure(2, ImportError("boom"))
    assert "pure-Python backend" not in str(record[0].message)


def test_backends_emptied_by_hand_still_explain_themselves(monkeypatch):
    """Tests simulate a missing extension by setting cppext._backends to
    None; explain() must still give a sentence then."""
    monkeypatch.setitem(cppext._backends, 3, None)
    monkeypatch.setitem(cppext._states, 3, dict(state="ok", error=None))
    assert "is not compiled" in cppext.explain(3)


def test_stale_extension_warns_once_loaded(fake_backend):
    d = _backend_folder(fake_backend)
    buildinfo.stamp(d)
    cppext._states[3] = dict(state="ok", error=None)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        cppext._warn_if_stale(3)  # sources unchanged: quiet
    assert cppext._states[3]["stale"] is False
    (fake_backend / "chain_session.h").write_text("// edited\n")
    with pytest.warns(RuntimeWarning, match="older than its C\\+\\+ sources"):
        cppext._warn_if_stale(3)
    assert cppext._states[3]["stale"] is True


def test_report_names_every_backend():
    text = doctor.report()
    for name in ("ITensor v3", "ITensor v2", "python", "julia_live",
                 "A chain with no itensor_version runs on"):
        assert name in text


@pytest.mark.parametrize("version", [2, 3])
def test_selftest_passes_on_a_compiled_backend(version):
    if not cppext.available(version):
        pytest.skip("ITensor v%s extension not compiled"%version)
    ok, msg = doctor.selftest(version)
    assert ok, msg


def test_selftest_of_a_missing_backend_fails_with_the_reason(monkeypatch):
    monkeypatch.setitem(cppext._backends, 2, None)
    ok, msg = doctor.selftest(2)
    assert not ok and "ITensor v2" in msg
