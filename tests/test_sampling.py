import numpy as np

from sctransform_cpu.sampling import (
    RMersenneTwister,
    SelectionManifest,
    r_nrd_density_weights,
)


def test_r45_uniform_stream_and_sample_are_exact():
    rng = RMersenneTwister(42)
    expected = np.asarray([
        0.91480604349635541,
        0.93707541329786181,
        0.28613953478634357,
        0.83044762606732547,
        0.64174551889300346,
        0.51909594913013279,
        0.73658831464126706,
        0.13466659723781049,
        0.65699229040183127,
        0.70506478403694928,
    ])
    actual = np.asarray([rng.uniform() for _ in range(len(expected))])
    assert np.array_equal(actual.view(np.uint64), expected.view(np.uint64))

    rng = RMersenneTwister(42)
    assert (rng.sample_uniform(240, 20) + 1).tolist() == [
        49, 229, 65, 153, 74, 228, 146, 122, 240, 128,
        47, 24, 71, 100, 89, 165, 110, 20, 154, 114,
    ]


def test_r45_weighted_density_sample_matches_r452():
    index = np.arange(1, 3001, dtype=np.float64)
    abundance = index * index / 1e6
    rng = RMersenneTwister(42)
    selected = rng.sample_weighted(r_nrd_density_weights(abundance), 2000) + 1
    assert selected[:30].tolist() == [
        899, 847, 2552, 1030, 1734, 2054, 1433, 2827, 1688, 1540,
        2197, 1492, 853, 2611, 2186, 839, 589, 2856, 2157, 1951,
        921, 2820, 626, 822, 2909, 2065, 2343, 918, 2221, 1,
    ]


def test_selection_manifest_round_trip(tmp_path):
    manifest = SelectionManifest(42, ("cell-a", "cell-b"), ("gene-z",))
    path = tmp_path / "selection.json"
    manifest.write(path)
    assert SelectionManifest.read(path) == manifest

