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
        
        print("LQM object:", lqm)
        print("Available methods:", [m for m in dir(lqm) if not m.startswith('_')])
        
        # Try getTool method
        try:
            # Let's see what tools are available
            print("\nTrying to get tools...")
            
            # Try different ways to get tools
            tool = lqm.getTool()
            print("Default tool:", tool)
            
        except Exception as e1:
            print("getTool() failed:", e1)
            
            # Try with parameters
            try:
                # Common tool names in BioLQM
                for tool_name in ['sbml', 'boolnet', 'ginml', 'petri']:
                    try:
                        tool = lqm.getTool(tool_name)
                        print(f"Tool {tool_name}:", tool)
                        break
                    except:
                        continue
                        
            except Exception as e2:
                print("getTool with params failed:", e2)
        
        # Try loadEngine
        try:
            engine = lqm.loadEngine()
            print("Engine:", engine)
            print("Engine methods:", [m for m in dir(engine) if not m.startswith('_')])
            
        except Exception as e3:
            print("loadEngine failed:", e3)
            
        # The key insight: let's try to access the service manager differently
        # Look for static methods or fields
        try:
            # Try accessing the class through JVM directly
            jvm = gateway.jvm
            
            # Check if we can get to BioLQM classes through different paths
            print("\nExploring BioLQM classes:")
            
            # Try different package paths
            for pkg_path in ['org.colomoto.biolqm', 'org.colomoto.biolqm.service']:
                try:
                    pkg = eval(f"jvm.{pkg_path}")
                    print(f"{pkg_path}: {pkg}")
                    print(f"  Type: {type(pkg)}")
                    if hasattr(pkg, '__dict__'):
                        print(f"  Contents: {dir(pkg)}")
                except:
                    continue
            
        except Exception as e4:
            print("Class exploration failed:", e4)
            
    except Exception as e:
        print("Overall error:", e)
        import traceback
        traceback.print_exc()
        
    proc.terminate()
else:
    print("Failed to get port")
    print("stderr:", proc.stderr.read())