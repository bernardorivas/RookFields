# test_entry_point.py
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
        # Use the entry point instead
        entry_point = gateway.entry_point
        print("Entry point:", entry_point)
        
        # Try the LQM method
        lqm = entry_point.LQM()
        print("LQM object:", lqm)
        print("LQM type:", type(lqm))
        print("LQM methods:", dir(lqm))
        
        # Try to access service manager through LQM
        if hasattr(lqm, 'service'):
            service = lqm.service()
            print("Service:", service)
            print("Service methods:", dir(service))
        
        # Try to get tools
        if hasattr(lqm, 'getTools'):
            tools = lqm.getTools()
            print("Tools:", tools)
        elif hasattr(lqm, 'tools'):
            tools = lqm.tools()
            print("Tools:", tools)
            
    except Exception as e:
        print("Error:", e)
        import traceback
        traceback.print_exc()
        
    proc.terminate()
else:
    print("Failed to get port")
    print("stderr:", proc.stderr.read())