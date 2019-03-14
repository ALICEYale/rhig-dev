/**
 * Attempt to reproduce the AliOADBContainer memory leak observed when running
 * the EMCal Correction Task. Unfortunately, this code is not sufficient to
 * reproduce the leak. Further information is available at this JIRA ticket:
 * https://alice.its.cern.ch/jira/projects/EMCAL/issues/EMCAL-135
 *
 * @author: Raymond Ehlers, Yale University <raymond.ehlers@cern.ch>
 */
#include <AliOADBContainer.h>
#include <AliDataFile.h>
#include <AliLog.h>

void aliOADBContainerLeakReproducer()
{
  AliLog::SetClassDebugLevel("AliDataFile", AliLog::kDebug);
  AliLog::SetClassDebugLevel("AliOADBContainer", AliLog::kDebug+2);
  AliOADBContainer * cont = new AliOADBContainer("");

  // Retrieve
  const std::string filename = AliDataFile::GetFileNameOADB("EMCAL/EMCALTemperatureCorrCalib.root");
  std::cout << "Filename: " << filename << "\n";
  cont->InitFromFile(filename.c_str(), "AliEMCALRunDepTempCalibCorrections");

  std::cout << "Sleep to catch up\n";
  sleep(3);

  std::cout << "Re-init\n";
  cont->InitFromFile(filename.c_str(), "AliEMCALRecalib");

  std::cout << "Deleting cont\n";
  delete cont;

}

int main()
{
  aliOADBContainerLeakReproducer();
  
  return 0;
}
