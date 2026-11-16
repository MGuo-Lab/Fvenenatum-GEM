import json
import subprocess

from fbaval.provenance import DEFAULT_REPOS, repo_state, write_provenance


def test_default_repos_is_this_repository():
    assert set(DEFAULT_REPOS) == {"Fvenenatum-GEM"}


def _make_repo(path):
    subprocess.run(["git", "init", "-q", str(path)], check=True)
    subprocess.run(["git", "-C", str(path), "-c", "user.name=t", "-c", "user.email=t@t",
                    "commit", "-q", "--allow-empty", "-m", "init"], check=True)


def test_repo_state_clean(tmp_path):
    _make_repo(tmp_path)
    state = repo_state(tmp_path)
    assert len(state["commit"]) == 40
    assert state["dirty"] is False


def test_repo_state_dirty(tmp_path):
    _make_repo(tmp_path)
    (tmp_path / "x.txt").write_text("x")
    assert repo_state(tmp_path)["dirty"] is True


def test_output_folders_do_not_mark_dirty(tmp_path):
    _make_repo(tmp_path)
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "x.csv").write_text("x")
    assert repo_state(tmp_path)["dirty"] is False


def test_repo_state_missing_repo(tmp_path):
    state = repo_state(tmp_path / "nope")
    assert state["commit"] is None


def test_write_provenance(tmp_path):
    repo = tmp_path / "repo"
    repo.mkdir()
    _make_repo(repo)
    out = tmp_path / "run1"
    path = write_provenance(out, repos={"gem": repo}, params={"o2_cap": 2.0})
    data = json.loads(path.read_text())
    assert data["repos"]["gem"]["dirty"] is False
    assert data["params"] == {"o2_cap": 2.0}
    assert "timestamp" in data
