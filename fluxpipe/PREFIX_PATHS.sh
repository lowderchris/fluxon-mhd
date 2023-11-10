#!/bin/zsh


# Define a function to check for the existence of CONFIG_FILE in specified directories
find_config_file() {
  local CONFIG_FILE="$1"
  # echo -e "\nLocating $CONFIG_FILE..."

  # Define an array of directories to search in
  local DIRECTORIES=("." ".." "fluxpipe")

  # Loop through the directories and check if the config.ini file exists
  for dir in "${DIRECTORIES[@]}"; do
    if [[ -f "$dir/$CONFIG_FILE" ]]; then
      # echo -e "\tFound $CONFIG_FILE in $dir/"
      # Output the path to the found file
      echo -e "$dir/$CONFIG_FILE"
      return  # Return without specifying a value
    fi
  done

  # If the file was not found in any of the directories, return an empty string
  return  # Return without specifying a value
}



CONFIG_FILE="config.ini"
# Call the function with the CONFIG_FILE argument
found_file_path=$(find_config_file $CONFIG_FILE)

# Check if the file was found and display the path if it was
if [[ -n "$found_file_path" ]]; then
  echo -e "The configuration file was found at: $found_file_path"
  # Continue with further actions here...
  # echo "\nUpdating paths using $found_file_path..."
else
  echo -e "The configuration file $CONFIG_FILE was not found in any of the specified directories."
  exit 1
fi

# Function to read and extract the value from config.ini based on the key passed
extract_path() {
  local key=$1
  local value=$(awk -F'=' -v key="$key" '$1 ~ key {print $2}' $found_file_path | sed 's/^\s*//;s/\s*$//')
  echo $(eval echo $value)  # Use eval to expand the tilde to $HOME if necessary
}

echo "Extracting paths..."
# Extract the paths from the config.ini file
SHELL_RC=$(extract_path "rc_path")
SHELL_RC_CUSTOM=$SHELL_RC"_custom"
RUN_SCRIPT=$(extract_path "run_script")
FL_PREFIX=$(extract_path "fl_prefix")
FL_MHDLIB=$(extract_path "fl_mhdlib")
DATA_DIR=$(extract_path "data_dir")
PYTHON_DIR=$(dirname $(which python3))
PL_PREFIX=$(dirname $(which perl))

echo Adjusting shell rc scripts...
# Create $SHELL_RC_CUSTOM file if it does not exist
if [[ ! -f "${SHELL_RC_CUSTOM}" ]]; then
  touch "${SHELL_RC_CUSTOM}"
  echo "Created ${SHELL_RC_CUSTOM} file."
fi

# Ensure that ~/.zshrc sources $SHELL_RC_CUSTOM
if ! grep -q "source ${SHELL_RC_CUSTOM}" "${SHELL_RC}"; then
  echo " " >> "${SHELL_RC}"
  echo "echo 'Loading ${SHELL_RC_CUSTOM}'" >> "${SHELL_RC}"
if ! grep -q "source ${SHELL_RC_CUSTOM}" "${SHELL_RC}"; then
  echo " " >> "${SHELL_RC}"
  echo "echo 'Loading ${SHELL_RC_CUSTOM}'" >> "${SHELL_RC}"
  # source ${SHELL_RC_CUSTOM}
  echo "source ${SHELL_RC_CUSTOM}" >> "${SHELL_RC}"
  echo "echo '\t\t Loaded ${SHELL_RC_CUSTOM}'" >> "${SHELL_RC}"
  echo "Added source command to ${SHELL_RC}"
  echo "source ${SHELL_RC_CUSTOM}" >> "${SHELL_RC}"
  echo "echo '\t\t Loaded ${SHELL_RC_CUSTOM}'" >> "${SHELL_RC}"
  echo "Added source command to ${SHELL_RC}"
else
  echo "\tThe RC File ${SHELL_RC} already sources ${SHELL_RC_CUSTOM}."
fi

# Function to update or set the environment variables
change_count=0
update_var() {
    local var_name=$1
    local new_value=$(eval echo $2)  # Use eval to expand the tilde to $HOME if necessary
    local current_value=$(grep "^export $var_name=" "$SHELL_RC_CUSTOM" | cut -d '=' -f2- | tr -d '"')

    # touch "$SHELL_RC_CUSTOM"
    # touch "$SHELL_RC"

    # Check and update the shell rc custom file
    if ! grep -q "^export $var_name=" "$SHELL_RC_CUSTOM"; then
        # If the line doesn't exist, append it
        echo "export $var_name=\"$new_value\"" >> "$SHELL_RC_CUSTOM"
        echo "Added $var_name to $SHELL_RC_CUSTOM:: $new_value"
        export $var_name=$new_value
        ((change_count+=1))
    elif [ "$current_value" = "$new_value" ]; then
        # If the value is the same, indicate no change
        a=1
        # echo "No change in $var_name in $SHELL_RC_CUSTOM:: $new_value"
    else
        # If the line exists but the value is different, replace it
        echo "$SHELL_RC_CUSTOM"
        # sed -i '' "s|^export $var_name=.*|export $var_name=\"$new_value\"|" "$HOME/.zshrc_custom"
        # sed -i '' "s|^export ${var_name}=.*|export ${var_name}=\"${new_value}\"|" "$HOME/.zshrc_custom"
        sed -i '' "s|^export ${var_name}=.*|export ${var_name}=\"${new_value}\"|" "$HOME/.zshrc_custom"

        echo "Updated $var_name in $SHELL_RC_CUSTOM to $new_value"
        export $var_name=$new_value
        ((change_count+=1))
    fi

    # Export the new value
    export $var_name="$new_value"
}

# Function to append a value to an environment variable if it isn't already present
append_var() {
    local var_name=$1
    local new_value=$(eval echo $2)  # Use eval to expand the tilde to $HOME if necessary
    local current_value=$(grep "^export $var_name=" "$SHELL_RC_CUSTOM" | cut -d '=' -f2- | tr -d '"')

    # Check if the value is already part of the variable
    if [[ ":$current_value:" != *":$new_value:"* ]]; then
        # If the value is not in current_value, append it
        new_value="${current_value:+$current_value:}$new_value"
        if ! grep -q "^export $var_name=" "$SHELL_RC_CUSTOM"; then
            # If the line doesn't exist, append it
            echo "export $var_name=\"$new_value\"" >> "$SHELL_RC_CUSTOM"
            # echo "Appended $var_name to $SHELL_RC_CUSTOM: $new_value"
        else
            # If the line exists, replace it with the new appended value
            sed -i '' "s|^export $var_name=.*|export $var_name=\"$new_value\"|" "$SHELL_RC_CUSTOM"
            # echo "Appended to $var_name in $SHELL_RC_CUSTOM: $new_value"
        fi
        # Export the new appended value
        export $var_name="$new_value"
    else
        # If the value is the same, indicate no change
        # echo "No change to $var_name in $SHELL_RC_CUSTOM: $new_value"
    fi
}

echo "Updating variables..."
echo $SHELL_RC_CUSTOM

# Update or set the extracted variables
update_var "SHELL_RC_CUSTOM" "$SHELL_RC_CUSTOM"
update_var "SHELL_RC" "$SHELL_RC"
update_var "RUN_SCRIPT" "$RUN_SCRIPT"
update_var "PYTHON_DIR" "$PYTHON_DIR"
update_var "PL_PREFIX" "$PL_PREFIX"
update_var "FL_PREFIX" "$FL_PREFIX"
update_var "FL_MHDLIB" "$FL_MHDLIB"
update_var "DATA_DIR" "$DATA_DIR"

echo "Updated $change_count environment variables."

add_subdirs_to_perl5lib_recursive() {
  local root_dir=$1  # Pass the base directory as an argument to the function

  if [[ -d "$root_dir" ]]; then  # Check if the directory exists
    while IFS= read -r -d '' dir; do
      append_var "PERL5LIB" $dir
      export PERL5LIB="$PERL5LIB:$dir"
    done < <(find "$root_dir" -type d -print0)

    echo "Subdirectories of $root_dir have been added to PERL5LIB."
  else
    echo "The directory $root_dir does not exist."
  fi
}

# Function to add immediate subdirectories of a given path to PERL5LIB, excluding certain directories
add_subdirs_to_perl5lib() {
  local root_dir=$1           # The base directory as an argument to the function
  shift                       # Shift the arguments so that $@ contains only the exclusion list
  local exclusions=("$@")     # Assign the rest of the arguments as an array of exclusions

  # Ensure the root directory exists
  if [[ ! -d "$root_dir" ]]; then
    echo "The directory $root_dir does not exist or is not accessible."
    return 1
  fi

  # Keep track of how many directories were added
  local added_count=0

  # Loop over all items within the root directory
  for dir in "$root_dir"/*; do
    # Check if the directory is in the exclusion list
    local exclude_dir=false
    for exclusion in "${exclusions[@]}"; do
      if [[ "$dir" == *"$exclusion"* ]]; then
        exclude_dir=true
        break
      fi
    done

    # If the item is a directory, not in exclusions, and not already in PERL5LIB, add it
    if [[ -d "$dir" && ":$PERL5LIB:" != *":$dir:"* && "$exclude_dir" == false ]]; then
      export PERL5LIB="$PERL5LIB:$dir"
      append_var "PERL5LIB" "$dir"
      ((added_count++))
      # echo "Added $dir to PERL5LIB"
    fi
  done

  if [[ $added_count -eq 0 ]]; then
    echo "No new subdirectories were added to PERL5LIB."
  else
    echo "Added $added_count new subdirectories to PERL5LIB."
  fi
}

export PERL5LIB=""
# add_subdirs_to_perl5lib "$FL_MHDLIB/fluxpipe/fluxpipe" "_Inline" "__pycache__"
export PERL5LIB="$PERL5LIB:$PL_PREFIX"
append_var "PERL5LIB" "$PL_PREFIX"

echo "Path Script Completed Successfully!\n"
