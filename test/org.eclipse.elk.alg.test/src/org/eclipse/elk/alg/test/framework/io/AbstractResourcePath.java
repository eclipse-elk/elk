/*******************************************************************************
 * Copyright (c) 2018, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.test.framework.io;

import java.io.File;
import java.io.FileFilter;
import java.util.ArrayList;
import java.util.List;

import com.google.common.base.Strings;

/**
 * Path to a resource to be used during tests. Resource paths can specify all files to be found in a directory (either
 * directly or recursively) or an actual file. They are based on a path specification relative to some base directory.
 * That base directory can be anything: a bundle's directory, a path specified through a system property, etc. The base
 * path is determined by the subclasses of this class. If a path specifies a whole number of files, they can be
 * filtered.
 * 
 * <p>
 * The relative path is interpreted relative to a base path. It can either refer to a concrete file, to all files in a
 * directory (in which case it must end in {@link #ALL_FILES_SUFFIX}), or to all files in a directory and its
 * sub directories (in which case it must end in {@link #ALL_FILES_RECURSIVE_SUFFIX}). Both the base path and the
 * relative path can simply use the slash {@code /} as their separator to be system-independent; those get converted
 * to the system-dependent path separator.
 * </p>
 * 
 * <p>
 * If a resource path describes the files in a directory (and possibly its sub directories), a list of all affected
 * resources can be obtained by calling 
 * </p>
 * 
 * <p>
 * In their constructor, subclasses must call {@link #initialize(String, String)} to initialize everything.
 * </p>
 */
public abstract class AbstractResourcePath {
    
    /** The suffix at the end of a relative path specification that targets all files in a directory. */
    public static final String ALL_FILES_SUFFIX = "/**";
    /**
     * The suffix at the end of a relative path specification that targets all files in a directory and its sub
     * directories.
     */
    public static final String ALL_FILES_RECURSIVE_SUFFIX = "/**/";

    /** The {@link File} that represents the actual resource. */
    private File file;
    /** Whether the resource path is a directory. */
    private boolean isDirectory = false;
    /** Whether files in the resource path should be looked up recursively. */
    private boolean isRecursive = false;
    /** The filter for looking up files in directories. */
    private FileFilter filter = null;
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Initialization
    
    /**
     * Sets up the path to this resource. The path is specified by a base path, which varies with each subclass, and a
     * relative path spec, which is resolved by this method relative to the base path.
     * 
     * @param base
     *            the base path, with or without trailing slash. Can be empty or {@code null}.
     * @param relative
     *            the relative path spec.
     */
    protected void initialize(final String base, final String relative) {
        if (relative == null) {
            throw new IllegalArgumentException("Relative path cannot be null.");
        }
        
        // Initialize base path
        String path = "";
        if (!Strings.isNullOrEmpty(base)) {
            path = base.replace("/", File.separator);
            if (!path.endsWith(File.separator)) {
                path += File.separator;
            }
        }
        
        // Relative path specification
        if (relative.endsWith(ALL_FILES_SUFFIX)) {
            // All files in the directory
            isDirectory = true;
            path += relative
                    .substring(0, relative.length() - ALL_FILES_SUFFIX.length())
                    .replace("/", File.separator);
            
        } else if (relative.endsWith(ALL_FILES_RECURSIVE_SUFFIX)) {
            // All files in the directory and its sub directories
            isDirectory = true;
            isRecursive = true;
            path += relative
                    .substring(0, relative.length() - ALL_FILES_RECURSIVE_SUFFIX.length())
                    .replace("/", File.separator);
            
        } else {
            // A single file
            path += relative.replace("/", File.separator);
        }
        
        file = new File(path);
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Accessors

    /**
     * Returns the {@link File} that corresponds to this resource path. The file may describe an actual file or a
     * directory.
     */
    public File getFile() {
        return file;
    }

    /**
     * Returns whether the specified path leads to a directory.
     */
    public boolean isDirectory() {
        return isDirectory;
    }

    /**
     * Returns whether to look for files recursively. This is only meaningful if {@link #isDirectory()} returns
     * {@code true}.
     */
    public boolean isRecursive() {
        return isRecursive;
    }

    /**
     * Returns the filter used to find files in directories.
     * 
     * @return the active file filter. May be {@code null}.
     */
    public FileFilter getFilter() {
        return filter;
    }
    
    /**
     * Sets the filter that should be used to look for files in directories. A {@code null} filter accepts all files.
     */
    public void setFilter(final FileFilter filter) {
        this.filter = filter;
    }

    /**
     * A method chaining version of {@link #setFilter(FileFilter)}.
     */
    public AbstractResourcePath withFilter(final FileFilter fileFilter) {
        setFilter(fileFilter);
        return this;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Resolver
    
    /**
     * Returns a list of resource paths that describe all of the actual resources described by this path. If this path
     * refers to a file, the list will contain only that file. If it refers to a directory, the list will contain all
     * of the files in the directory that are accepted by the filter. If it is set to be recursive, sub directories are
     * traversed as well.
     */
    public List<AbsoluteResourcePath> listResources() {
        List<AbsoluteResourcePath> result = new ArrayList<>();
        
        if (file.exists()) {
            if (isDirectory && file.isDirectory()) {
                // We're a directory, and we're also supposed to be one. Let another method gather all files
                fillFileList(file, result);
                
            } else if (!isDirectory && file.isFile()) {
                // We're not a directory, just a humble file, so add us to the result
                result.add(new AbsoluteResourcePath(file.getAbsolutePath()));
            }
        }
        
        return result;
    }
    
    private void fillFileList(final File directory, final List<AbsoluteResourcePath> result) {
        // List children (null filter accepts everything)
        File[] children = directory.listFiles(filter);
        
        for (File child : children) {
            if (child.isFile()) {
                // Files are always added to the results
                result.add(new AbsoluteResourcePath(child.getAbsolutePath()));
                
            } else if (child.isDirectory() && isRecursive) {
                // Recurse into the sub directory
                fillFileList(child, result);
            }
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Stuff
    
    @Override
    public String toString() {
        return file.getAbsolutePath();
    }
    
}
