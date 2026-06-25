#-----------------------------------------------------------------------------
set(MODULE_NAME PelvicFloor)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  ${MODULE_NAME}Lib/__init__
  ${MODULE_NAME}Lib/GlobalUtils
  ${MODULE_NAME}Lib/PelvicFloorUtils
  ${MODULE_NAME}Lib/SlicerScripts
  ${MODULE_NAME}Lib/UserScripts
  ${MODULE_NAME}Lib/MeshRegistration
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
