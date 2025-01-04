options("repos" = c(CRAN="https://mirrors.aliyun.com/CRAN/"))

# install.packages("nadiv")
# install.packages("optparse")
# install.packages("AlphaSimR")


# 获取已安装包的信息
installed_packages <- installed.packages()

# 获取指定包的版本
package_version <- installed_packages["nadiv", "Version"]
# 打印指定包的版本
print(package_version)

# 获取指定包的版本
package_version <- installed_packages["AlphaSimR", "Version"]
# 打印指定包的版本
print(package_version)

# 获取指定包的版本
package_version <- installed_packages["optparse", "Version"]
# 打印指定包的版本
print(package_version)
