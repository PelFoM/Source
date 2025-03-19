#-----------------------------------------------------------------------------
set(MODULE_NAME PelvicFloor)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  ${MODULE_NAME}Lib/__init__
  ${MODULE_NAME}Lib/PelvicFloorUtil
  ${MODULE_NAME}Lib/ModelTools
  ${MODULE_NAME}Lib/MeshMorphing
  ${MODULE_NAME}Lib/UserScripts
  ${MODULE_NAME}Lib/Utils
  #TemplateModel/__init__
  #TemplateModel/Landmarks
  )

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_RESOURCES
  #Resources/Icons/SlicerAdvancedGear-Small.png
  #${MODULE_NAME}Lib/Resources/UI/settings.ui
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )
