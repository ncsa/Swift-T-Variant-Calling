package provide myextension 0.0

namespace eval myextension {
  proc double { x } {
      return [ expr $x * 2 ]
  }
}
