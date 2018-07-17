// check if this task should be performed

bool check_task(auto name) {
  auto input = InputGroup("tasks.in","tasks");
  auto dotask = input.getYesNo(name,false);
  return dotask ;
}

auto get_int_value(auto name) {
  auto input = InputGroup("tasks.in","tasks");
  auto N = input.getInt(name,1);
  return N ;
}




auto get_float_value(auto name) {
  auto input = InputGroup("tasks.in","tasks");
  auto N = input.getReal(name,1);
  return N ;
}




auto get_str(auto name) {
  auto input = InputGroup("tasks.in","tasks");
  return input.getString(name,"");
}



