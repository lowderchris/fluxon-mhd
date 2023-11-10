#!/bin/zsh

# Ensure the script stops on first error
set -e

[[ ! $(pwd) =~ fluxpipe ]] && cd fluxpipe

echo "\n\nInstalling FluxPipe..."

sudo apt-get update
sudo apt install zsh gcc gnuplot-qt make


# # Source the .zshrc_custom
# if [[ -f $SHELL_RC_CUSTOM ]]; then
#   source $SHELL_RC_CUSTOM
# fi

# Define the desired Perl version
desired_perl_version="perl-5.32.0"

# Check if perlbrew is installed
if ! command -v perlbrew &>/dev/null; then
    echo "\tInstalling perlbrew..."
    \curl -L https://install.perlbrew.pl | bash
    # source ~/perl5/perlbrew/etc/bashrc

    # The line to add to the shell configuration file
    PERLBREW_INIT_LINE='source ~/perl5/perlbrew/etc/bashrc'

    # $PERLBREW_INIT_LINE

    # Check if the line is already in the configuration file
    if ! grep -qxF "$PERLBREW_INIT_LINE" "$SHELL_RC_CUSTOM"; then
        # If the line isn't present, add it to the shell configuration file
        echo "$PERLBREW_INIT_LINE" >> "$SHELL_RC_CUSTOM"
        echo "Added perlbrew initialization line to $SHELL_RC_CUSTOM"
    else
        echo "Perlbrew initialization line already present in $SHELL_RC_CUSTOM"
    fi

fi

# Check if the desired Perl version is available in perlbrew
if perlbrew list | grep -q "$desired_perl_version"; then
    echo "\t$desired_perl_version is already installed."
else
    echo "\tInstalling $desired_perl_version with perlbrew..."
    perlbrew install "$desired_perl_version"
fi

# Always switch to the desired Perl version
perlbrew switch "$desired_perl_version"

# Install cpanminus if not already present
if ! command -v cpanm &>/dev/null; then
    echo "\tInstalling cpanminus..."
    yes | perlbrew install-cpanm
fi

# Print which perl is currently active
echo "Current Perl: $(which perl)"

# # Check if Homebrew is installed
# if ! command -v brew &>/dev/null; then
#     /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
#     sudo brew install gnuplot
#     sudo brew install fftw
#     sudo brew install qt

# else
#     echo "\tHomebrew already installed!"
# fi

# Check if conda command exists
if command -v conda &>/dev/null; then
    # Check if conda is initialized for zsh
    if [ ! -f $SHELL_RC ] || ! grep -q "conda.sh" $SHELL_RC; then
        conda init zsh
        conda init bash
        source $SHELL_RC
        # You can add conda-specific operations here if needed
    else
        echo "\tAnaconda already initialized in $SHELL_RC!"
    fi

    # Check if the fluxenv conda environment exists
    conda_envs=$(conda info --envs)
    if ! echo "$conda_envs" | grep -q "fluxenv"; then
        echo "\nPython Environment Creation:"
        # Your code to create the fluxenv if it doesn't exist
        yes | conda create --name fluxenv --file requirements-conda.txt
        # yes | conda create -n fluxenv
        conda activate fluxenv
        yes | conda install pip
        yes | pip install -e .
        yes | pip install -r requirements-pip.txt
    fi

else
    echo "\tConda is not installed. Falling back to virtualenv and pip."

    # Check if virtualenv is installed, and if not, install it
    if ! command -v virtualenv &>/dev/null; then
        yes | sudo apt update
        yes | sudo apt upgrade
        yes | sudo apt install python3-pip
        yes | sudo apt install python3-virtualenv
        # pip install virtualenv
    fi

# Define the desired virtual environment path
venv_path="$FL_PREFIX/fluxenv_pip"

# Check if the virtual environment already exists
if [ ! -d "$venv_path" ]; then
    # Create a new virtual environment if it doesn't exist
    virtualenv "$venv_path"

    # Activate the virtual environment
    source "$venv_path/bin/activate"

    # Install packages from the two requirements.txt files
    pip install -r requirements-pip.txt
    # Install packages from requirements-conda.txt and continue even if it fails
    pip install -r requirements-conda.txt || true
    # Install the current project in the virtual environment
    "$venv_path/bin/pip" install -e .
fi
    # Activate the existing virtual environment
    source "$venv_path/bin/activate"


fi

source $SHELL_RC

add_subdirs_to_paths_recursive() {
  local root_dir=$1  # Pass the base directory as an argument to the function
    local depth=${2:-3}

  if [[ -d "$root_dir" ]]; then  # Check if the directory exists
    while IFS= read -r -d '' dir; do
      append_var "PERL5LIB" $dir
      export PERL5LIB="$PERL5LIB:$dir"
      export PATH="$PATH:$dir"
    done < <(find "$root_dir" -maxdepth $depth -type d -print0)

    echo "Subdirectories of $root_dir have been added to PERL5LIB."
  else
    echo "The directory $root_dir does not exist."
  fi
}

# Function to append a value to an environment variable if it isn't already present
append_var() {
    local var_name="$1"
    local new_value="$(eval echo "$2")"  # Use eval to expand the tilde to $HOME if necessary

    # Extract the current value using awk
    local current_value=$(awk -F'=' -v var_name="$var_name" '$1 == "export " var_name {gsub(/"/, "", $2); print $2}' "$SHELL_RC_CUSTOM")

    # Check if the value is already part of the variable
    if [[ ":$current_value:" != *":$new_value:"* ]]; then
        # If the value is not in current_value, append it
        new_value_cat="${current_value:+$current_value:}$new_value"
        if ! grep -q "^export $var_name=" "$SHELL_RC_CUSTOM"; then
            # If the line doesn't exist, append it
            echo "export $var_name=\"$new_value_cat\"" >> "$SHELL_RC_CUSTOM"
            echo "Appended $var_name to $SHELL_RC_CUSTOM: $new_value"
        else
            # If the line exists, replace it with the new appended value
            sed -i'' "s|^export $var_name=.*|export $var_name=\"$new_value_cat\"|" "$SHELL_RC_CUSTOM"
            echo "Appended to $var_name in $SHELL_RC_CUSTOM: $new_value"
        fi
        # Export the new appended value
        export "$var_name"="$new_value"
    else
        # If the value is the same, indicate no change
        # echo "No change to $var_name in $SHELL_RC_CUSTOM: $new_value"
    fi
}


change_count=0
update_var() {
    local var_name=$1
    local new_value=$(eval echo $2)  # Use eval to expand the tilde to $HOME if necessary
    local current_value=$(grep "^export $var_name=" "$SHELL_RC_CUSTOM" | cut -d '=' -f2- | tr -d '"')
    # Check and update the shell rc custom file
    if ! grep -q "^export $var_name=" "$SHELL_RC_CUSTOM"; then
        # If the line doesn't exist, append it
        echo "export $var_name=\"$new_value\"" >> "$SHELL_RC_CUSTOM"
        echo "Added $var_name to $SHELL_RC_CUSTOM:: $new_value"
        ((change_count+=1))
    elif [ "$current_value" = "$new_value" ]; then
        # If the value is the same, indicate no change
        a=1
        # echo "No change in $var_name in $SHELL_RC_CUSTOM:: $new_value"
    else
        # If the line exists but the value is different, replace it
        sed -i '' "s|^export $var_name=.*|export $var_name=\"$new_value\"|" "$SHELL_RC_CUSTOM"
        echo "Updated $var_name in $SHELL_RC_CUSTOM to $new_value"
        ((change_count+=1))
    fi

    # Export the new value
    export $var_name="$new_value"
}

add_subdirs_to_paths_recursive $FL_PREFIX
add_subdirs_to_paths_recursive $PL_PREFIX 4

# conda activate fluxenv

# make it so zshenv calls zshrc
# grep -qF "source ~/.zshrc" ~/.zshenv || echo "source ~/.zshrc" >> ~/.zshenv

# echo the current directory


# Check and install PDL and other dependencies
cd $PL_PREFIX
if ! perl -MPDL -e 1 &>/dev/null; then
    cpanm PDL
fi
# export PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX/lib/perl5"
export PERL_MM_OPT="INSTALL_BASE=$PL_PREFIX"
sudo cpan local::lib
if ! grep -Fxq 'eval "$(perl -I$PL_PREFIX/lib/perl5 -Mlocal::lib=$PL_PREFIX)"' $SHELL_RC_CUSTOM; then
        echo '\neval "$(perl -I$PL_PREFIX/lib/perl5 -Mlocal::lib=$PL_PREFIX)"' >> $SHELL_RC_CUSTOM
fi
sudo cpan File::ShareDir
sudo cpan PDL::Graphics::Gnuplot
sudo cpan Math::RungeKutta
# ... Add other dependencies as needed ...

# Install FLUXpipe pdl
# (
# cpanm --installdeps .
# cpanm --notest --local-lib=$FL_PREFIX .
# if ! grep -q 'export PERL5LIB=$FL_PREFIX/lib/perl5:$PERL5LIB' $SHELL_RC_CUSTOM; then
#     echo 'export PERL5LIB=$FL_PREFIX/lib/perl5:$PERL5LIB' >> $SHELL_RC_CUSTOM
# fi
# )

cd $FL_MHDLIB/fluxpipe
rm -rf blib/

cd $FL_MHDLIB
ln -sf fluxpipe/fluxpipe/runners/config_runner.py config_runner.py
echo "\n\nFluxpipe installation complete!\n\n"

if command -v conda &>/dev/null; then
    echo "Remember to run 'conda activate fluxenv' before running 'python fluxpipe/fluxpipe/runners/config_runner.py'\n\n"
else
    echo "Remember to run 'source fluxenv_pip/bin/activate' before running 'python fluxpipe/fluxpipe/runners/config_runner.py'\n\n"
fi

echo -n "Do you want to test the config_runner? ([yes]/no): "
read response
if [[ "$response" == "yes" || "$response" == "y" || "$response" == "" ]]; then
    echo "Running the file..."
    if command -v conda &>/dev/null; then
        conda activate fluxenv
    else
        # Activate the virtual environment
        source $FL_PREFIX/fluxenv_pip/bin/activate
    fi
    echo "Environment Activated"
    cd $FL_MHDLIB
    python config_runner.py
fi
