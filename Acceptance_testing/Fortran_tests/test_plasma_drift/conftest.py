# conftest.py file to setup pytest
import pytest

# Command line arguments
def pytest_addoption(parser):
    parser.addoption("--trun", action="store", default="0",help="--trun: 0/1")
    parser.addoption("--ttest", action="store", default="0",help="--ttest: 0/1")
    parser.addoption("--tpath", action="store", default='',help="--ttest: 0/1")
    
# Pass arguments to the test function
@pytest.fixture
def trun(request):
    return request.config.getoption("--trun") 
@pytest.fixture
def ttest(request):
    return request.config.getoption("--ttest",default='1')   
@pytest.fixture
def tpath(request):
    return request.config.getoption("--tpath",default='') 
