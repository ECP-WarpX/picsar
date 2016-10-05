# conftest.py file to setup pytest
import pytest

# Command line arguments
def pytest_addoption(parser):
    parser.addoption("--trun", action="store", default="1",help="--trun: 0/1")
    parser.addoption("--ttest", action="store", default="1",help="--ttest: 0/1")
    parser.addoption("--tshow", action="store", default="1",help="--tshow: 0/1")
    parser.addoption("--tpath", action="store", default='',help="--tpath: path where the test is run")
    
# Pass arguments to the test function
@pytest.fixture
def trun(request):
    return request.config.getoption("--trun",default='1') 
@pytest.fixture
def ttest(request):
    return request.config.getoption("--ttest",default='1')   
@pytest.fixture
def tpath(request):
    return request.config.getoption("--tpath",default='') 
@pytest.fixture
def tshow(request):
    return request.config.getoption("--tshow",default='1') 
