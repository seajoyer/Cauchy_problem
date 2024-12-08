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
          pname = "ode_solver";
          version = "0.1.0";
          name = "ode_solver";

          src = ./cpp;

          nativeBuildInputs = with pkgs; [ gnumake ];

          buildInputs = with pkgs; [ gnuplot ];

          buildPhase = ''
            g++ -std=c++17 -lm main.cpp -o ode_solver
          '';

          installPhase = ''
            mkdir -p $out/bin
            cp ode_solver $out/bin/
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
              ${cppProject}/bin/ode_solver
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
