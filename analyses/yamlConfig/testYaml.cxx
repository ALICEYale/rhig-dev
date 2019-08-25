#include <string>
#include <iostream>

#include <AliLog.h>
#include <AliYAMLConfiguration.h>

bool CompareStringValues(const std::string & inputString, const std::string & retrievedString)
{
  bool returnValue = false;
  if (inputString != retrievedString) {
    // Admittedly, this is both redundant and useless, since it will go to the fatal next
    returnValue = false;
    AliFatalGeneralF("AliYAMLConfigurationTest", "Input string \"%s\" does not match retrieved string \"%s\"", inputString.c_str(), retrievedString.c_str());
  }
  else {
    std::cout << "Success!\n";
    returnValue = true;
  }

  return returnValue;
}

void testYaml()
{
  // Just the SetClassDebugLevel was sufficient to work with AliLog
  //AliLog::SetGlobalLogLevel(AliLog::kDebug);
  //AliLog::SetGlobalDebugLevel(5);
  AliLog::SetClassDebugLevel("AliYAMLConfiguration", 5);
  //AliLog::SetModuleDebugLevel("AliYAMLConfiguration", 5);
  AliDebugGeneral("AliYAMLConfiguration", 0, "Test");
  AliDebugGeneralStream("AliYAMLConfiguration", 2) << "Test stream\n";
  AliDebugGeneralStream(__PRETTY_FUNCTION__, 2) << "Test pretty function\n";
  PWG::Tools::AliYAMLConfiguration testConfig;
  std::string retrievedValue = "";
  std::string testStr = "testStr";
  std::string testStr2 = "testStr2";

  // Named defined in the "test.yaml" file
  std::string configName = "hi";
  testConfig.AddConfiguration("test.yaml");
  testConfig.AddConfiguration("test2.yaml", "name2");
  testConfig.Print(AliInfoGeneralStream("AliYAMLConfigurationTest"));

  std::cout << "Basic reading tests\n";
  retrievedValue = "";
  testConfig.GetProperty("name", retrievedValue);
  CompareStringValues("hi", retrievedValue);

  std::cout << "Reading overritten value test\n";
  retrievedValue = "";
  testConfig.GetProperty("hello2:world2", retrievedValue);
  CompareStringValues(testStr2, retrievedValue);

  std::cout << "Initial writing test\n";
  testConfig.WriteProperty("hello:world", testStr, configName);
  // Check the value
  retrievedValue = "";
  testConfig.GetProperty("hello:world", retrievedValue);
  CompareStringValues(testStr, retrievedValue);

  std::cout << "Test existing YAML Node\n";
  testConfig.WriteProperty("hello2:world", testStr2, configName);
  // Check the value
  retrievedValue = "";
  testConfig.GetProperty("hello2:world", retrievedValue);
  CompareStringValues(testStr2, retrievedValue);

  std::cout << "Test vector\n";
  std::vector<std::string> test = { "Hi", "There", "World" };
  testConfig.WriteProperty("hello2:testVector", test);
  // Check the value
  std::vector<std::string> retrievedVector;
  testConfig.GetProperty("hello2:testVector", retrievedVector);
  for (unsigned int i = 0; i < test.size(); i++)
  {
    if (test.at(i) != retrievedVector.at(i)) {
      AliFatalGeneralF("AliYAMLConfigurationTest", "Vector values mismatch! index: %u, input: %s, retrieved: %s", i, test.at(i).c_str(), retrievedVector.at(i).c_str());
    }
  }
  // If we get to this point, then we've succeeded.
  std::cout << "Success!\n";

  std::cout << "Test complex object\n";
  std::map<std::string, std::vector<std::string>> testMap;
  testMap["testMapKey"] = test;
  testConfig.WriteProperty("hello3", testMap);
  // Check the value
  std::map<std::string, std::vector<std::string>> retrievedMap;
  testConfig.GetProperty("hello3", retrievedMap);
  for (auto pair : testMap)
  {
    std::vector<std::string> tempVector = retrievedMap.at(pair.first);
    for (unsigned int i = 0; i < pair.second.size(); i++)
    {
      if (pair.second.at(i) != tempVector.at(i)) {
        AliFatalGeneralF("AliYAMLConfigurationTest", "Vector values in map mismatch! index: %u, input: %s, retrieved: %s", i, pair.second.at(i).c_str(), tempVector.at(i).c_str());
      }
    }
  }
  // If we get to this point, then we've succeeded.
  std::cout << "Success!\n";

  // Test editing a node directly
  std::cout << "Testing adding a value directly to the node\n";
  auto configPair = testConfig.GetConfiguration(configName);
  std::string directConfigStr = "directConfigValue";
  configPair.second[directConfigStr] = directConfigStr;

  // Print for debug information
  testConfig.Print(AliInfoGeneralStream("AliYAMLConfigurationTest"), configName);

  // Check the value
  retrievedValue = "";
  testConfig.GetProperty(directConfigStr, retrievedValue);
  CompareStringValues(directConfigStr, retrievedValue);

  std::cout << "\n*******************************************\n";
  std::cout << "* Success!\n";
  std::cout << "*******************************************\n\n";
  std::cout << "Final configurations\n";
  std::cout << testConfig;
}
