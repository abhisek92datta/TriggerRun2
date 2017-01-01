//Add includes for your classes here
#include "DataFormats/Common/interface/Wrapper.h"
#include "TriggerRun2/TriggerAnalyzer/interface/TriggerStudyEventVars.h"

#include <vector>

namespace {
   struct TriggerRun2_TriggerAnalyzer{
  
     //add 'dummy' Wrapper variable for each class type you put into the Event
     triggerStudyEventVars triggerStudyEventVarsDummy0;
     std::vector<triggerStudyEventVars> triggerStudyEventVarsDummy1;  
     edm::Wrapper<triggerStudyEventVars> triggerStudyEventVarsDumm2;
     edm::Wrapper<std::vector<triggerStudyEventVars> > triggerStudyEventVarsDummy4;
   };
}
