{
  description = "Cauchy problem";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = nixpkgs.legacyPackages.${system};

        pythonEnv = pkgs.python3.withPackages
          (ps: with ps; [ numpy matplotlib sympy tqdm ]);

        cppProject = pkgs.stdenv.mkDerivation {
          pname = "cauchy_problem_example";
          version = "0.1.0";
          name = "cauchy_problem-0.1.0";

          src = ./cpp;

          nativeBuildInputs = with pkgs; [ cmake gnumake ];

          buildInputs = with pkgs; [ gnuplot ];

          cmakeFlags = [
            "-DCMAKE_BUILD_TYPE=Release"
            "-DCMAKE_CXX_FLAGS=-std=c++17"
          ];

          buildPhase = ''
            cmake .
            make VERBOSE=1 -j $NIX_BUILD_CORES
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp cauchy_problem_example $out/bin/
          '';

          fixupPhase = ''
            patchelf --set-rpath "${pkgs.lib.makeLibraryPath [ pkgs.gnuplot ]}:${pkgs.stdenv.cc.cc.lib}/lib" $out/bin/cauchy_problem_example
          '';
        };

        pythonProject = pkgs.stdenv.mkDerivation {
          pname = "python-cauchy_problem";
          version = "0.1.0";
          name = "python-cauchy_problem-0.1.0";

          src = ./py;

          nativeBuildInputs = [ pythonEnv ];

          installPhase = ''
            mkdir -p $out/bin $out/lib/python
            cp -r . $out/lib/python/
            cat > $out/bin/run-python <<EOF
            #!${pythonEnv}/bin/python
            import sys
            import os

            sys.path.insert(0, '$out/lib/python')
            sys.path.insert(0, '$out/lib/python/cauchy_problem')

            import demo
            demo.main()
            EOF
            chmod +x $out/bin/run-python
          '';
        };

      in {
        packages = {
          cpp = cppProject;
          py = pythonProject;
          default = cppProject;
        };

        apps = {
          cpp = flake-utils.lib.mkApp {
            drv = pkgs.writeShellScriptBin "run-cpp" ''
              export PATH="${pkgs.gnuplot}/bin:$PATH"
              ${cppProject}/bin/cauchy_problem_example
            '';
          };
          py = flake-utils.lib.mkApp {
            drv = pythonProject;
            name = "run-python";
          };
          default = self.apps.${system}.cpp;
        };

        devShells.default = pkgs.mkShell {
          name = "dev_shell";

          nativeBuildInputs = with pkgs; [
            cmake
            cmake-language-server
            bear
            gnumake
            clang
            libcxx
            gnuplot
            ccache
            git
            pyright
            pythonEnv
          ];

          buildInputs = with pkgs; [ boost catch2 ];

          shellHook = ''
            export CXXFLAGS="''${CXXFLAGS:-} -I${pkgs.catch2}/include -std=c++17"

            export CCACHE_DIR=$HOME/.ccache
            export PATH="$HOME/.ccache/bin:$PATH"

            alias c=clear

            echo "======================================"
            echo "$(cmake   --version | head -n 1)"
            echo "$(g++     --version | head -n 1)"
            echo "$(make    --version | head -n 1)"
            echo "$(python  --version | head -n 1)"
            echo "$(gnuplot --version | head -n 1)"
            echo ""
            echo "Build the project:  nix build"
            echo "Run C++ project:    nix run .#cpp"
            echo "Run Python project: nix run .#py"
            echo ""
            echo "Happy coding!"
          '';
        };
      });
}
