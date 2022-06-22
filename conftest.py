from armicontrib.armiopenmc.tests import openmcTestingApp


def pytest_sessionstart(session):
    import armi

    if not armi.isConfigured():
        armi.configure(openmcTestingApp.OpenmcTestingApp())
