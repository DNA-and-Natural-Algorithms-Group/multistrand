
function reset() { block = 0; py_block = 0; flush = 0; buf0 = ""; buf1 = "" }
BEGIN { reset() }

! /^###/    # keep non-comment lines
  /^###/ {  # un-comment commented blocks referring to PyObject methods
    line = substr($0, 4)
    switch (line) {
        case "{": { block++; break }
        case "}": { if (block) flush++; break }
        default:  { if (block && index(line, "fun:_PyObject_")) py_block++ }
    }
    if (block)  { buf0 = buf0 $0; buf1 = buf1 line }
    if (!flush) { if (block) { buf0 = buf0 ORS; buf1 = buf1 ORS } else print }
    else        { if (py_block) { print buf1 } else { print buf0 }; reset()  }
}
