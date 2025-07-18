import subprocess
import time
from py4j.java_gateway import JavaGateway, GatewayParameters

# Start GINsim
print("Starting GINsim...")
proc = subprocess.Popen(['GINsim', '-py'], 
                       stdout=subprocess.PIPE, 
                       stderr=subprocess.PIPE,
                       text=True)

time.sleep(2)
port_line = proc.stdout.readline().strip()
print(f"Port: {port_line}")

if port_line.isdigit():
    port = int(port_line)
    gateway = JavaGateway(gateway_parameters=GatewayParameters(port=port))
    
    print("Gateway connected!")
    
    try:
        # Use the entry point
        entry_point = gateway.entry_point
        lqm = entry_point.LQM()
        
        # Test load methods
        print("\nTesting load methods:")
        
        # Create a simple test model file
        test_sbml = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1" 
      xmlns:qual="http://www.sbml.org/sbml/level3/version1/qual/version1" qual:required="true">
  <model>
    <qual:listOfQualitativeSpecies>
      <qual:qualitativeSpecies qual:id="A" qual:maxLevel="1"/>
      <qual:qualitativeSpecies qual:id="B" qual:maxLevel="1"/>
    </qual:listOfQualitativeSpecies>
  </model>
</sbml>"""
        
        with open('test_model.sbml', 'w') as f:
            f.write(test_sbml)
        
        # Try to load the model
        try:
            model = lqm.load('test_model.sbml')
            print("Model loaded:", model)
            print("Model type:", type(model))
            if model:
                print("Model methods:", [m for m in dir(model) if not m.startswith('_')])
        except Exception as e:
            print("Load failed:", e)
            
        # Try loadModel
        try:
            model = lqm.loadModel('test_model.sbml')
            print("loadModel result:", model)
        except Exception as e:
            print("loadModel failed:", e)
        
        # Test formats
        print("\nTesting formats:")
        try:
            format_obj = lqm.getFormat('sbml')
            print("SBML format:", format_obj)
            if format_obj:
                print("Format methods:", [m for m in dir(format_obj) if not m.startswith('_')])
        except Exception as e:
            print("getFormat failed:", e)
            
        # Test tools  
        print("\nTesting tools:")
        tool_names = ['stable', 'trace', 'trapspace', 'export']
        for tool_name in tool_names:
            try:
                tool = lqm.getTool(tool_name)
                print(f"Tool '{tool_name}':", tool)
                if tool:
                    print(f"  Methods:", [m for m in dir(tool) if not m.startswith('_')])
            except Exception as e:
                print(f"getTool('{tool_name}') failed:", e)
                
        # Test modifier
        print("\nTesting modifiers:")
        try:
            modifier = lqm.getModifier('reduction')
            print("Reduction modifier:", modifier)
        except Exception as e:
            print("getModifier failed:", e)
            
        # Try to access service manager directly through JVM
        print("\nTrying direct JVM access:")
        try:
            jvm = gateway.jvm
            # Try to get LQMServiceManager
            service_manager = jvm.org.colomoto.biolqm.service.LQMServiceManager
            print("ServiceManager:", service_manager)
            
            # Try getServices
            try:
                services = service_manager.getServices()
                print("Services:", services)
            except:
                pass
                
            # Try getTools
            try:
                tools = service_manager.getTools()
                print("Tools:", tools)
            except:
                pass
                
        except Exception as e:
            print("Direct JVM access failed:", e)
            
    except Exception as e:
        print("Overall error:", e)
        import traceback
        traceback.print_exc()
        
    finally:
        proc.terminate()
        import os
        if os.path.exists('test_model.sbml'):
            os.remove('test_model.sbml')
else:
    print("Failed to get port")
    print("stderr:", proc.stderr.read())