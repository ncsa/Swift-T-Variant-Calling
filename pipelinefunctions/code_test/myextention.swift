// Run this file using:
//     swift-t -r $PWD myextention.swift

@pure
(int o) double (int i) "myextension" "0.0" [
  "set <<o>> [ myextension::double <<i>> ]"
];


trace("\n\n" + 3 +" \t " + double(3)+ "\n\n");
