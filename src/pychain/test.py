# this script checks that everythong is fine, and
# that all the libraries are present


try:
  import entanglement

except:
  print("entanglement library not properly compiled")
  exit()
