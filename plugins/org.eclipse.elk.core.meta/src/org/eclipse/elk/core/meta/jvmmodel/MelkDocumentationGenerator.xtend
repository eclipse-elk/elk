/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.jvmmodel

import com.google.common.base.CharMatcher
import com.google.common.collect.Iterables
import java.io.PrintWriter
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.util.EnumSet
import java.util.LinkedList
import org.eclipse.elk.core.meta.MetaDataRuntimeModule.MelkOutputConfigurationProvider
import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdBundle
import org.eclipse.elk.core.meta.metaData.MdBundleMember
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdModel
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.elk.graph.properties.GraphFeature
import org.eclipse.xtext.common.types.JvmEnumerationType
import org.eclipse.xtext.common.types.JvmField
import org.eclipse.xtext.common.types.JvmParameterizedTypeReference
import org.eclipse.xtext.common.types.JvmTypeReference
import org.eclipse.xtext.generator.AbstractFileSystemAccess2
import org.eclipse.xtext.generator.IFileSystemAccess
import org.eclipse.xtext.nodemodel.util.NodeModelUtils
import org.eclipse.xtext.xbase.XExpression
import org.eclipse.xtext.xbase.compiler.JvmModelGenerator
import org.eclipse.elk.core.meta.metaData.MdOptionDependency

/**
 * Generates documentation from the source model that pertains to a *.melk file.
 * <p>
 * For each algorithm, option and group in the model, a Markdown file is created in the directory specified by
 * the vm argument {@code elk.metadata.documentation.outputPath}. This path usually points to the docs folder in the
 * git repository for elk, when the Oomph setup appends the vm argument to the eclipse.ini.
 * The docs folder in the master branch will be hosted on a Github page.
 * </p>
 * 
 * @author dag
 */
class MelkDocumentationGenerator extends JvmModelGenerator {
    /**
     * The {@code IFileSystemAccess} used to read additional documentation files stored within the Eclipse project that
     * contains the model.
     */
    IFileSystemAccess fsa
    /** The place where the generated algorithm documentation is stored. */
    Path algorithmsOutputPath
    /** The place where the generated layout option documentation is stored. */
    Path optionsOutputPath
    /** The place where the generated layout option group documentation is stored. */
    Path optionGroupsOutputPath
    /** The place where images are stored. */
    Path imageOutputPath
    /** The directory containing additional documentation files. */
    Path projectDocumentationSourceFolder
    
    /**
     * The method {@code internalDoGenerate} is called for each {@link MdModel} derived from a *.melk file.
     * 
     * @param model
     *      the {@link MdModel} that documentation is generated for.
     * @param fsa
     *      the {@code IFileSystemAccess} used to read files
     */
    def dispatch void internalDoGenerate(MdModel model, IFileSystemAccess fsa) {
        this.fsa = fsa;
        
        // The property is set as a vm argument in the eclipse.ini and holds the output folder for documentation.
        val propertyOutputPath = System.getProperty("elk.metadata.documentation.outputPath")
        // If the property is not set, don't generate any documentation.
        if (propertyOutputPath === null) {
            return
        }

        // Make sure the output paths exist and are empty
        setupOutputPaths(propertyOutputPath);
        
        // If the model has no name or no bundle section, nothing can be generated.
        if (model.name === null || model.bundle === null) {
            return
        }
        
        val bundle = model.bundle
        val members = bundle.members
        
        // The {@code documentationFolder} in the bundle section tells us where to look for additional 
        // documentation files.
        if (bundle.documentationFolder !== null) {
            projectDocumentationSourceFolder = Paths.get(bundle.documentationFolder)
        } else {
            // The default folder for additional documentation is a docs folder on project level.
            projectDocumentationSourceFolder = Paths.get("docs")
        }

        // write documentation for layout algorithms, groups and options
        for (algorithm : members.filter(MdAlgorithm)) {
            algorithm.writeDoc(algorithmsOutputPath)
        }
        
        for (option : members.allOptionDefinitions) {
            option.writeDoc(optionsOutputPath)
        }
        
        for (group : members.allGroupDefinitions) {
            group.writeDoc(optionGroupsOutputPath)
        }
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Input / Output
    
    /**
     * Sets up the documentation output paths. Normally, we would try and purge the directories first to be sure that
     * there is no old content. However, in a typical automatic build this generator will be invoked multiple times. If
     * each invocation purged the output directories, there wouldn't be much documentation left.
     * 
     * @param docsOutputFolder the folder the documentation should be generated into.
     */
    private def setupOutputPaths(String docsOutputPathString) {
        val Path docsOutputPath = Paths.get(docsOutputPathString);
        
        // Inside the documentation folder, we need to obtain access to the content/reference subfolder
        val Path referencePath = docsOutputPath.resolve("content").resolve("reference");
        
        // Algorithms folder (remove and create to purge the thing)
        algorithmsOutputPath = referencePath.resolve("algorithms");
        Files.createDirectories(algorithmsOutputPath);
        
        // Layout options folder (remove and create to purge the thing)
        optionsOutputPath = referencePath.resolve("options");
        Files.createDirectories(optionsOutputPath);
        
        // Layout option groups folder (remove and create to purge the thing)
        optionGroupsOutputPath = referencePath.resolve("groups");
        Files.createDirectories(optionGroupsOutputPath);
        
        // Images folder (remove and create to purge the thing)
        imageOutputPath = docsOutputPath.resolve("static").resolve("img_gen");
        Files.createDirectories(imageOutputPath);
    }
    
    /** 
     * This is a wrapper for {@link #writeDoc(String, String) writeDoc(String, String)} 
     * that uses the full name of a member as filename.
     * 
     * @param member
     *      the {@link MdBundleMember} that documentation is generated for
     * @param outputPath
     *      the folder the bundle member's documentation is to be generated into
     */
    private def void writeDoc(MdBundleMember member, Path outputPath) {
        val fileName = member.qualifiedName.toHugoIdentifier
        writeDoc(fileName, outputPath, member.generateDoc)
    }
    
    /**
     * Creates a file in {@code outputPath} to store the generated Markdown in.
     * 
     * @param fileName
     *      name of the file
     * @param outputPath
     *      the folder the bundle member's documentation is to be generated into
     * @param documentation
     *      generated markdown string
     */
    private def void writeDoc(String fileName, Path outputPath, String documentation) {
        var file = fileName
        // add file extension for Markdown
        if (!file.endsWith(".md")) {
            file += ".md"
        }

        var PrintWriter out
        try {
            out = new PrintWriter(outputPath.resolve(file).toString);
            out.print(documentation)
        } catch (Exception exception) {
            exception.printStackTrace
        } finally {
            out?.close
        }
    }
    
    /**
     * Copies images or any file from within the project to the place where the generated documentation is located.
     * 
     * @param path
     *      the path to the image relative to the eclipse project
     * @param newFilename
     *      file name for the copy (should be associable to the issuing member)
     */
    private def copyImageToOutputPath(String path, String newFileName) {
        if (path === null || newFileName === null) {
            return
        }

        try {
            var iStream = (fsa as AbstractFileSystemAccess2).readBinaryFile(path,
                MelkOutputConfigurationProvider.AD_INPUT)
            Files.copy(iStream, imageOutputPath.resolve(newFileName), StandardCopyOption.REPLACE_EXISTING)
        } catch (Exception exception) {
            exception.printStackTrace
        }
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Algorithm Documentation

    /**
     * The dispatch method {@code generateDoc} is called for each {@link MdBundleMember} that is supposed to be
     * documented.
     * This generates a documentation String in Markdown format.
     * 
     * @param algorithm
     *      the {@link MdAlgorithm} from the {@link MdModel}
     * @return
     *      the Markdown String
     */
    private def dispatch String generateDoc(MdAlgorithm algorithm) {
        // The frontmatter header (between dashes) tells the static-site-generator how to render this document.
        var doc = '''
        ---
        title: "«algorithm.label ?: algorithm.name»"
        displayid: "«algorithm.qualifiedName»"
        menu:
          main:
            identifier: "alg-«algorithm.qualifiedName.toHugoIdentifier»"
            parent: "Algorithms"
        ---
        
        '''
        
        if (algorithm.previewImage !== null) {
            val newFileName = algorithm.qualifiedName.toHugoIdentifier + "_preview_" 
                                + algorithm.previewImage.substring(algorithm.previewImage.lastIndexOf('/') + 1)
                                
            // copy previewImage into images folder within the rest of the documentation
            algorithm.previewImage.copyImageToOutputPath(newFileName)
            doc += "{{< image src=\"" + newFileName + "\" alt=\"Preview Image\" gen=\"1\" >}}\n\n"
        }
        
        doc += '''
        Property | Value
        -------- | -----
        *Identifier:* | `«algorithm.qualifiedName»`
        *Meta Data Provider:* | `«algorithm.bundle.targetClass»`
        
        '''
        
        if (!algorithm.description.isNullOrEmpty) {
            doc += '''
            ## Description
            
            «algorithm.description.trimNewlineTabsAndReduceToSingleSpace»
            
            '''
        }
        
        // If the algorithm belongs to a category, get its description as well.
        if (algorithm.category !== null) {
            doc += '''
            ## Category: «algorithm.category.label ?: algorithm.category.name»
            
            '''
            
            if (!algorithm.category.description.isNullOrEmpty) {
                doc += '''
                «algorithm.category.description.trimNewlineTabsAndReduceToSingleSpace»
                
                '''
            }
        }
    
        // table of supported graph features
        if (!algorithm.supportedFeatures.empty) {
            doc += '''
            ## Supported Graph Features
            
            Name | Description
            ----|----
            '''
            for (supportedFeature : algorithm.supportedFeatures) {
                doc += supportedFeature.literal.split("_").map[it.toFirstUpper].join(" ") + 
                        // the description of each feature is retrieved from core.options.GraphFeature
                        " | " + GraphFeature.valueOf(supportedFeature.name()).description + "\n"
            }
            doc += "\n"
        }
        
        if (algorithm.documentation !== null) {
            doc += '''
            ## Additional Documentation
            
            «algorithm.documentation.additionalDocumentation(algorithm.qualifiedName.toHugoIdentifier)»
            
            '''
        }
        
        // table of supported options
        if (!algorithm.supportedOptions.empty) {
            doc += '''
            ## Supported Options
            
            Option | Default Value
            ----|----
            '''
            for (supportedOption : algorithm.supportedOptions.sortBy[(it.option.label ?: it.option.name).toLowerCase]) {
                var optionFileName = supportedOption.option.qualifiedName.toHugoIdentifier
                
                // If the supported option specifies additional documentation, it is specific to this algorithm,
                // therefore an alternative file is created for the option documentation that differs from the general
                // documentation.
                if (supportedOption.documentation !== null) {
                    optionFileName = supportedOption.option.qualifiedName.toHugoIdentifier + "_" +
                        algorithm.qualifiedName.toHugoIdentifier
                    writeDoc(
                        optionFileName,
                        optionsOutputPath,
                        supportedOption.option.generateDoc(
                            (supportedOption.option.label ?: supportedOption.option.name) + " ("
                                + (algorithm.label ?: algorithm.name) + ")",
                            optionFileName,
                            supportedOption.documentation.additionalDocumentation(optionFileName))
                    )
                }
                
                val optionFileLink = '''{{< relref "reference/options/«optionFileName».md" >}}''';
                
                doc += '''[«supportedOption.option.label ?: supportedOption.option.name»](«optionFileLink»)''' +
//                ''' | `«supportedOption.option.type?.simpleName/*?.replace("<","&lt;").replace(">","&gt;")*/»`''' +
                ''' | `«supportedOption.value?.text ?: supportedOption.option.defaultValue.text»`''' +
//                ''' | «supportedOption.option.qualifiedName.replace(".", "&#8203;.")»''' +
                "\n"
            }
            doc += "\n"
        }
        
        return doc
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Option Documentation
    
    /**
     * The dispatch method {@code generateDoc} is called for each {@link MdBundleMember} that is supposed to be
     * documented.
     * This generates a documentation String in Markdown format.
     * 
     * @param option
     *      the {@link MdOption} from the {@link MdModel}
     * @return
     *      the Markdown String
     */
    private def dispatch String generateDoc(MdOption option) {
        generateDoc(option, option.label ?: option.name,
            option.qualifiedName.toHugoIdentifier,
            option.documentation.additionalDocumentation(option.qualifiedName.toHugoIdentifier))
    }
    
    /**
     * This method is wrapped by {@link #generateDoc(MdOption)}.
     * Here the title and additional documentation can be specified so that an algorithm that supports this option can 
     * create an alternate version of this documentation, that is specific to the algorithm.
     * 
     * @param option
     *      the {@link MdOption} from the {@link MdModel}
     * @param title
     *      the title for the documentation page
     * @param id
     *      the identifier of the options page to be generated
     * @param additionalDoc
     *      a String of additional documentation that will be appended 
     * @return
     *      the Markdown String
     * @see #additionalDocumentation
     */
    private def String generateDoc(MdOption option, String title, String id, String additionalDoc) {
        // The frontmatter header (between dashes) tells the static-site-generator how to render this document.
        var doc = '''
        ---
        title: "«title»"
        displayid: "«option.qualifiedName»"
        menu:
          main:
            identifier: "option-«id»"
            parent: "LayoutOptions"
        ---
        
        '''
        
        val kinds = #[
            "deprecated"->option.deprecated,
            "advanced"->option.advanced,
            "programmatic"->option.programmatic,
            "output"->option.output,
            "global"->option.global].filter[it.value]
        
        doc += '''
        
        Property | Value
        -------- | -----
        «if (!kinds.empty) '''*Type:* | «kinds.map[it.key].join(", ")»'''»
        *Identifier:* | `«option.qualifiedName»`
        *Meta Data Provider:* | `«option.bundle.targetClass»`
        *Value Type:* | `«option.type?.identifier»`«if (option.type?.type instanceof JvmEnumerationType) " (Enum)"»
        «if (option.type.possibleValues !== null) "*Possible Values:* | " + 
            option.type.possibleValues.map[
                "`" + it.simpleName + "`"
                + it.annotations?.sortBy[it.annotation?.simpleName]
                        .join(" (", " and ", ")",["*`@" + it.annotation?.simpleName + "`*"])
                ].join("<br>")»
        «if (option.defaultValue !== null) '''*Default Value:* | `«option.defaultValue.text
                    »` (as defined in `«option.bundle.idPrefix ?: (option.bundle.eContainer as MdModel).name»`)'''»
        «if (option.lowerBound !== null) '''*Lower Bound:* | `«option.lowerBound.text»`'''»
        «if (option.upperBound !== null) '''*Upper Bound:* | `«option.upperBound.text»`'''»
        «if (!option.targets.empty) "*Applies To:* | " + option.targets.map[it.literal].join(", ")»
        «if (!option.legacyIds.empty) "*Legacy Id:* | " + option.legacyIds.map(["`" + it + "`"]).join(", ")»
        «if (!option.dependencies.empty) "*Dependencies:* | " + option.dependencies.map[toOptionDocsString].join(", ")»
        «if (!option.groups.empty) "*Containing Group:* | " + option.groups.map[toOptionDocsString].join(" -> ")»
        «if (option.description !== null) "\n### Description\n\n" + option.description.trimNewlineTabsAndReduceToSingleSpace»
        «if (!additionalDoc.nullOrEmpty) "\n## Additional Documentation\n\n" + additionalDoc»
        '''

        return doc
    }
    
    /**
     * Generates a proper String representation for the given option dependency.
     */
    private def String toOptionDocsString(MdOptionDependency dependency) {
        // Link to the layout option we depend upon
        var result = "[" + dependency.target.qualifiedName + "]({{< relref \"reference/options/"
                + dependency.target.qualifiedName.toHugoIdentifier + ".md\" >}})";
        
        if (dependency.value !== null) {
            result += " (`" + getText(dependency.value) + "`)";
        }
        
        return result;
    }
    
    /**
     * Generates a proper String representation for the given group.
     */
    private def toOptionDocsString(MdGroup group) {
        return "[" + group.name+ "]({{< relref \"reference/groups/" + group.qualifiedName.toHugoIdentifier + ".md\" >}})"
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Group Documentation
    
    /**
     * The dispatch method {@code generateDoc} is called for each {@link MdBundleMember} that is supposed to be
     * documented.
     * This generates a documentation String in Markdown format.
     * 
     * @param group
     *      the {@link MdGroup} from the {@link MdModel}
     * @return
     *      the Markdown String
     */
    private def dispatch String generateDoc(MdGroup group) {
        val title = new StringBuffer(group.groups.map[it.name].join("."));
        if (title.length == 0) {
            title.append(group.name);
        } else {
            title.append(".").append(group.name);
        }
        
        var doc = '''
        ---
        title: "«title.toString()»"
        displayid: "«group.qualifiedName»"
        menu:
          main:
            identifier: "group-«group.qualifiedName.toHugoIdentifier»"
            parent: "LayoutOptionGroups"
        ---
        
        Property | Value
        -------- | -----
        *Identifier:* | `«group.qualifiedName»`
        
        «if (!group.children.filter(MdOption).empty) "\n## Options\n\n" + group.children.filter(MdOption).map[
            "* [" + (it.label ?: it.name) + "]({{< relref \"reference/options/" + it.qualifiedName.toHugoIdentifier + ".md\" >}})"].join("\n")»
        «if (!group.children.filter(MdGroup).empty) "\n## Subgroups\n\n" + group.children.filter(MdGroup).map[
            "* [" + it.name + "]({{< relref \"reference/groups/" + it.qualifiedName.toHugoIdentifier + ".md\" >}})"].join("\n")»
        «if (group.documentation !== null) "\n## Additional Documentation\n\n" 
            + group.documentation.additionalDocumentation(group.qualifiedName.toHugoIdentifier)»
        '''
        return doc
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Additional Documentation
    
    /**
     * Extracts additional documentation from file or text.
     * A leading {@code @} signals that the following string is supposed to be interpreted as a file name.
     * The file has to be located at the path specified by the {@code documentationFolder} in the bundle section, that
     * defaults to "docs".
     * <p>
     * Images that are linked to in Markdown syntax are copied to the {@code outputPath} and their URLs are adjusted
     * accordingly.
     * </p>
     * 
     * @param documentation
     *      {@code @<file name>} or text in Markdown format
     * @param fileNamePrefix
     *      prefix of the file name for the copy to associate it with the member it belongs to
     * @return
     *      the Markdown String with replaced image references
     */
    private def String additionalDocumentation(String documentation, String fileNamePrefix) {
        if (documentation === null) {
            return ""
        }
        
        var doc = ""
        
        if (!documentation.startsWith("@")) {
            doc = documentation.trimNewlineTabsAndReduceToSingleSpace
        } else {
            try {
                // AbstractFileSystemAccess2 needs an OutputConfiguration to work, that is defined in MetaDataRuntimeModule
                doc = (fsa as AbstractFileSystemAccess2).readTextFile(
                    projectDocumentationSourceFolder.resolve(documentation.substring(1)).toString,
                    MelkOutputConfigurationProvider.AD_INPUT).toString
            } catch (Exception exception) {
                exception.printStackTrace
            }
        }
        
        // copy all images to the output location and update their URLs in the documentation
        // regex doesn't work here because links may contain parentheses that have to be counted
        var res = ""
        var int i = 0
        var int l = doc.length
        while (i < l) {
            res += doc.charAt(i)
            i++
            // ![ indicates the beginning of an image reference
            if (res.endsWith("![")) {
                // Skip over alt text between braces. We are interested in the path surrounded by parentheses.
                while (i < l && !res.endsWith("](")) {
                    res += doc.charAt(i)
                    i++
                }
                res = res.substring(0, res.lastIndexOf("!["));
                var path = ""
                var int openedParens = 1
                // count parentheses to determine the closing one
                while (i < l && openedParens > 0) {
                    path += doc.charAt(i)
                    i++
                    switch (doc.charAt(i).toString) { // xtend generates a String from ')' not a char
                    	case "(": openedParens++
                    	case ")": openedParens--
                    }
                }
                // remove whitespace and possible title statement from path
                var imgTitle = ""
                path = path.trim
                if (path.endsWith('\"')) {
                    imgTitle = " " + path.substring(path.lastIndexOf('\"', path.length - 2))
                    path = path.substring(0, path.lastIndexOf('\"', path.length - 2)).trim
                }
                // copy file into images folder within the rest of the documentation
                val newFileName = fileNamePrefix + "_" + path.substring(path.lastIndexOf('/') + 1)
                projectDocumentationSourceFolder.resolve(path).toString.copyImageToOutputPath(newFileName)
                // replace the URL with the new path
                res += "{{< image src=\"" + newFileName + "\" alt=\"" + imgTitle + "\" gen=\"1\" >}}\n\n"
                i++;
            }
        }
        
        return res
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // String Processing
    
    /**
     * The qualified name of a member consists of the bundles idPrefix, the containing groups and the members name
     * joined by dots.
     */
    private def String getQualifiedName(MdBundleMember member) {
        val bundle = member.bundle
        val model = bundle.eContainer as MdModel
        var prefix = bundle.idPrefix ?: model.name

        if (prefix.endsWith(member.name)) {
           prefix = prefix.substring(0, prefix.lastIndexOf('.')) 
        }
        
        return prefix
               + (if (member.groups.empty) '' else '.')
               + member.groups.map[it.name].join('.') 
               + '.' + member.name
    }
    
    /**
     * Turns the given option identifier into an identifier to be used by Hugo. This includes replacing dots by
     * hyphens, and ensuring that no identifier ends with "index" since that causes Hugo to do strange things.
     * 
     * @param optionName the option name to turn into an identifier fit for use with Hugo.
     * @return the Hugo-compatible identifier.
     */
    private def String toHugoIdentifier(String optionName) {
        if (optionName.endsWith(".index")) {
            // If the last component of the name ends in "index", Hugo doesn't generate links properly (#268)
            return (optionName + "2").replace('.', '-');
        } else {
            return optionName.replace('.', '-');
        }
    }
    
    /**
     * Readable text representation of an XExpression.
     */
    private def String getText(XExpression exp) {
        if (exp === null) {
            return "<not defined>"
        }
        var text = NodeModelUtils.getNode(exp).text.trim
        // If it's a function name, append parentheses. We just infer that from the lowercaseness of the first letter
        // from the last segment, if there are dot separated segments at all.
        val fName = text.substring(text.lastIndexOf('.') + 1)
        if (Character.isLowerCase(fName.charAt(0)) && fName != text && !fName.endsWith(")")) {
            text += "()"
        }
        return text
    }
    
    /**
     * Turns the string into something we can safely include in Markdown.
     */
    private def String trimNewlineTabsAndReduceToSingleSpace(String string) {
        CharMatcher.breakingWhitespace().replaceFrom(string, ' ').replaceAll(" +", " ")
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Object Tree Traversal
    
    /**
     * Possible values for enumeration types.
     */
    private def Iterable<? extends JvmField> getPossibleValues(JvmTypeReference type) {
        if (type !== null) {
            val jvmType = type.type
            if (jvmType instanceof JvmEnumerationType) {
                val enumType = jvmType
                return enumType.literals
            // in case of EnumSet the inner type is extracted
            } else if (jvmType.identifier == EnumSet.canonicalName) {
                val enumType = (type as JvmParameterizedTypeReference).arguments.head.type as JvmEnumerationType
                return enumType.literals
            }
        }
        return null
    }
    
    /**
     * Iterator for nested groups (flattened).
     */
    private def Iterable<MdGroup> getAllGroupDefinitions(Iterable<? extends MdBundleMember> elements) {
        val groups = elements.filter(MdGroup)
        return Iterables.concat(groups, groups.map[it.children.getAllGroupDefinitions].flatten)
    }
    
    /**
     * Iterator for nested options (flattened).
     */
    private def Iterable<MdOption> getAllOptionDefinitions(Iterable<? extends MdBundleMember> elements) {
        return Iterables.concat(elements.filter(MdOption),
                                elements.filter(MdGroup).map[it.children.getAllOptionDefinitions].flatten)
    }
    
    /**
     * Retrieve bundle associated with a member.
     */
    private def MdBundle getBundle(MdBundleMember member) {
        var parent = member.eContainer
        while (!(parent instanceof MdBundle)) {
            parent = parent.eContainer
        }
        return parent as MdBundle
    }
    
    /**
     * Groups containing this member.
     */
    private def Iterable<MdGroup> getGroups(MdBundleMember member) {
        val groups = new LinkedList
        var group = member.eContainer
        while (group instanceof MdGroup) {
            groups.addFirst(group) // ordered from outer to inner
            group = group.eContainer
        }
        return groups
    }
}
