# -*- Python -*-
# Drone CI Starlark configuration file.
# https://docs.drone.io/pipeline/scripting/starlark/
# Run `drone starlark convert --stdout` to verify `.drone.star`.
def get_test_commands():
  return [
    "mkdir -p build",
    "cd build",
    'CFLAGS="$ARCH_FLAGS $DIAGNOSTIC_FLAGS" CXXFLAGS="$ARCH_FLAGS $DIAGNOSTIC_FLAGS" meson .. || (cat meson-logs/meson-log.txt; false)',
    '"$(command -v ninja || command -v ninja-build)" -v test'
  ]

def get_apt_install_commands(extra_pkgs = []):
  return [
    "apt-get -yq update",
    "apt-get -yq install %s ninja-build git-core python3-pip parallel" % " ".join(extra_pkgs),
    "pip3 install meson",
  ]

def get_dnf_install_commands(extra_pkgs = []):
  return [
    "dnf install -y %s ninja-build git-core python3-pip parallel findutils" % " ".join(extra_pkgs),
    "pip3 install meson",
  ]

def get_yum_install_commands(extra_pkgs = []):
  return [
    "yum install -y epel-release",
    "yum install -y %s meson ninja-build git-core parallel" % " ".join(extra_pkgs),
  ]

def get_default_job():
  return {
    "kind": "pipeline",
    "type": "docker",
    "name": "",
    "platform": {
      "os": "linux",
    },
    "steps": [
      {
        "name": "test",
        "image": "ubuntu:bionic",
        "environment": {},
        "commands": [
        ]
      }
    ],
    "custom": {
      "before_install": [
        "uname -m",
        "cat /proc/cpuinfo",
      ],
      "install": [],
      "before_script": [
        "git submodule --quiet update --init --recursive",
      ],
      "script": get_test_commands()
    },
    "trigger": {
      "branch": {
        "exclude": [
          "master"
        ]
      }
    }
  }

def get_jobs():
  job_clang9_armv7 = {
    "name": "clang-9 armv7",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "clang-9",
          "CXX": "clang++-9",
          "ARCH_FLAGS": "-march=armv7a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["clang-9"])
    }
  }

  job_clang9_armv8 = {
    "name": "clang-9 armv8",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "clang-9",
          "CXX": "clang++-9",
          "ARCH_FLAGS": "-march=armv8a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["clang-9"])
    }
  }

  job_gcc8_armv7 = {
    "name": "gcc-8 armv7",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "gcc-8",
          "CXX": "g++-8",
          "ARCH_FLAGS": "-march=armv7-a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Wextra -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["gcc-8", "g++-8"])
    }
  }

  job_gcc8_armv8 = {
    "name": "gcc-8 armv8",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "gcc-8",
          "CXX": "g++-8",
          "ARCH_FLAGS": "-march=armv8-a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Wextra -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["gcc-8", "g++-8"])
    }
  }

  job_clang7_armv7 = {
    "name": "clang-7 armv7",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "clang-7",
          "CXX": "clang++-7",
          "ARCH_FLAGS": "-march=armv7a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["clang-7"])
    }
  }

  job_clang7_armv8 = {
    "name": "clang-7 armv8",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "clang-7",
          "CXX": "clang++-7",
          "ARCH_FLAGS": "-march=armv8a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["clang-7"])
    }
  }

  job_gcc7_armv7 = {
    "name": "gcc-7 armv7",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "gcc-7",
          "CXX": "g++-7",
          "ARCH_FLAGS": "-march=armv7-a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Wextra -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["gcc-7", "g++-7"])
    }
  }

  job_gcc7_armv8 = {
    "name": "gcc-7 armv8",
    "platform": {
      "arch": "arm",
    },
    "steps": [
      {
        "environment": {
          "CC": "gcc-7",
          "CXX": "g++-7",
          "ARCH_FLAGS": "-march=armv8-a -mfpu=neon",
          "DIAGNOSTIC_FLAGS": "-Wextra -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["gcc-7", "g++-7"])
    }
  }

  job_clang9_aarch64 = {
    "name": "clang-9 aarch64",
    "platform": {
      "arch": "arm64",
    },
    "steps": [
      {
        "environment": {
          "CC": "clang-9",
          "CXX": "clang++-9",
          "ARCH_FLAGS": "-march=armv8a+simd+crypto+crc",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_apt_install_commands(["clang-9"])
    }
  }

  job_fedora = {
    "name": "fedora",
    "steps": [
      {
        "image": "fedora:latest",
        "environment": {
          "CC": "gcc",
          "CXX": "g++",
          "ARCH_FLAGS": "-march=native",
          "DIAGNOSTIC_FLAGS": "-Wextra -Werror",
        }
      }
    ],
    "custom": {
      "install": get_dnf_install_commands(["gcc", "gcc-c++"])
    }
  }

  job_fedora_clang_arm64_flags = {
    "name": "fedora clang arm64 flags",
    "platform": {
      "arch": "arm64",
    },
    "steps": [
      {
        "image": "fedora:latest",
        "environment": {
          "CC": "clang",
          "CXX": "clang++",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        }
      }
    ],
    "custom": {
      "install": get_dnf_install_commands(["clang"]),
      "script": [
        # optflags RPM macro works with gcc.
        # Some flags and specs are not available with clang.
        # https://lists.fedoraproject.org/archives/list/packaging@lists.fedoraproject.org/message/W5UFLUADNB4VF3OBUBSNAPOQL6XBCP74/
        "ARCH_FLAGS=$(rpm -E '%{optflags}' | sed -e 's| -fstack-clash-protection||' -e 's| -specs=[^ ]*||g')",
      ] + get_test_commands()
    }
  }

  job_centos7_clang3 = {
    "name": "centos7 clang3",
    "steps": [
      {
        "image": "centos:7",
        "environment": {
          "CC": "clang",
          "CXX": "clang++",
          "DIAGNOSTIC_FLAGS": "-Weverything -Werror",
        },
        "failure": "ignore"
      }
    ],
    "custom": {
      # gcc, gcc-c++ are necessary to build on clang.
      "install": get_yum_install_commands(["clang", "gcc", "gcc-c++"])
    }
  }

  return [
    job_clang9_armv7,
    job_clang9_armv8,
    job_gcc8_armv7,
    job_gcc8_armv8,
    job_clang7_armv7,
    job_clang7_armv8,
    job_gcc7_armv7,
    job_gcc7_armv8,
    job_clang9_aarch64,
    # job_fedora,
    # job_fedora_clang_arm64_flags,
    job_centos7_clang3,
  ]

def main(ctx):
  merged_jobs = []
  for job in get_jobs():
    out = get_default_job()

    # Merge the each elements in the dict.
    for key, value in job.items():
      if type(value) == "list":
        for index, item in enumerate(value):
          out[key][index].update(item)
      elif type(value) == "dict":
        out[key].update(value)
      else:
        out[key] = value

    # Create commands list from custom elements.
    for element in ["before_install", "install", "before_script", "script"]:
      out["steps"][0]["commands"].extend(out["custom"][element])

    # Remove unused custom element.
    out.pop("custom", None)

    merged_jobs.append(out)

  return merged_jobs
