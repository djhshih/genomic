with import <nixpkgs> {};
stdenv.mkDerivation rec {
	name = "env";
	env = buildEnv { name = name; paths = buildInputs; };
	buildInputs = [ cmake lua boost mesa SDL freetype ftgl ];
}
