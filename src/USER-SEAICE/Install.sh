# Install/unInstall USER-SEAICE package files
# mode = 0/1/2 for uninstall/install/update

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test ! -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files

for file in *.cpp *.h; do
    action $file
done
