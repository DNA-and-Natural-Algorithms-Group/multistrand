
# Documentation: https://devenv.sh/

{ pkgs, lib, config, inputs, ... }:

let
  buildInputs = with pkgs; [ stdenv.cc.cc libuv zlib ];
in

{
  env = {
    LD_LIBRARY_PATH = "${with pkgs; lib.makeLibraryPath buildInputs}";
  };

  packages = with pkgs; [
    git
  ];

  languages.python = {
    enable = true;
    package = pkgs.python313;
    uv = {
      enable = true;
      sync.enable = true;
    };
  };

  enterShell = ''
    . .devenv/state/venv/bin/activate
  '';
}
