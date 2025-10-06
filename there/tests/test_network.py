import pytest
import there.network as tn

@pytest.fixture
def sample_network():
    # Create a sample network for testing
    return {
        'nodes': [1, 2, 3],
        'edges': [(1, 2), (2, 3)],
        'weights': [5, 10]
    }

def test_prepare_cycle_net(sample_network):

    # Test the prepare_cycle_net function
    cycle_net = tn.prepare_cycle_net(sample_network)
    assert 'cycle' in cycle_net
    assert len(cycle_net['nodes']) == len(sample_network['nodes'])
    assert len(cycle_net['edges']) == len(sample_network['edges'])
