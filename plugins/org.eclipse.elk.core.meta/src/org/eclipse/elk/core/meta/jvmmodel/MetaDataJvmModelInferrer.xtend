/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.jvmmodel

import com.google.common.collect.Iterables
import com.google.inject.Inject
import java.util.EnumSet
import java.util.LinkedList
import java.util.List
import org.eclipse.elk.core.meta.metaData.MdAlgorithm
import org.eclipse.elk.core.meta.metaData.MdBundle
import org.eclipse.elk.core.meta.metaData.MdBundleMember
import org.eclipse.elk.core.meta.metaData.MdCategory
import org.eclipse.elk.core.meta.metaData.MdGroup
import org.eclipse.elk.core.meta.metaData.MdModel
import org.eclipse.elk.core.meta.metaData.MdOption
import org.eclipse.elk.core.meta.metaData.MdOptionDependency
import org.eclipse.elk.core.meta.metaData.MdOptionSupport
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.properties.Property
import org.eclipse.xtend2.lib.StringConcatenationClient
import org.eclipse.xtext.common.types.JvmDeclaredType
import org.eclipse.xtext.common.types.JvmEnumerationType
import org.eclipse.xtext.common.types.JvmGenericType
import org.eclipse.xtext.common.types.JvmParameterizedTypeReference
import org.eclipse.xtext.common.types.JvmPrimitiveType
import org.eclipse.xtext.common.types.JvmTypeReference
import org.eclipse.xtext.common.types.JvmVisibility
import org.eclipse.xtext.common.types.util.Primitives
import org.eclipse.xtext.util.Strings
import org.eclipse.xtext.xbase.jvmmodel.AbstractModelInferrer
import org.eclipse.xtext.xbase.jvmmodel.IJvmDeclaredTypeAcceptor
import org.eclipse.xtext.xbase.jvmmodel.JvmTypesBuilder

/**
 * Infers a JVM model from the source model. 
 *
 * <p>The JVM model should contain all elements that would appear in the Java code which is generated from
 * the source model. Other models link against the JVM model rather than the source model.</p>     
 */
class MetaDataJvmModelInferrer extends AbstractModelInferrer {

    /** Convenience API to build and initialize JVM types and their members. */
    @Inject extension JvmTypesBuilder
    /** Convenience API to wrap primitives. */
    @Inject extension Primitives
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // MAIN METADATA PROVIDER CLASS GENERATION

    /**
     * The dispatch method {@code infer} is called for each instance of the
     * given element's type that is contained in a resource.
     * 
     * @param model
     *            the model to create one or more {@link JvmDeclaredType declared types} from.
     * @param acceptor
     *            each created {@link JvmDeclaredType type} without a container should be passed to the acceptor
     *            in order to get attached to the current resource. The acceptor's
     *            {@link IJvmDeclaredTypeAcceptor#accept(org.eclipse.xtext.common.types.JvmDeclaredType) accept(..)}
     *            method takes the constructed empty type for the pre-indexing phase. This one is further
     *            initialized in the indexing phase using passed closure.
     * @param isPreIndexingPhase
     *            whether the method is called in a pre-indexing phase, i.e.
     *            when the global index is not yet fully updated. You must not
     *            rely on linking using the index if isPreIndexingPhase is
     *            <code>true</code>.
     */
    def dispatch void infer(MdModel model, IJvmDeclaredTypeAcceptor acceptor, boolean isPreIndexingPhase) {
        // If we don't have a model name or the bundle section, don't bother
        if (model.name === null || model.bundle === null) {
            return
        }
        val bundle = model.bundle
        
        // Create a class for the general option definitions to be used by clients as well as for the registration
        // method that instantiates the ILayoutMetaDataProviders
        acceptor.accept(bundle.toClass(bundle.qualifiedTargetClass)) [
            superTypes += typeRef("org.eclipse.elk.core.data.ILayoutMetaDataProvider")
            fileHeader = model.documentation
            documentation = bundle.documentation
            
            // 1. Public constants for all declared layout options
            for (property : bundle.members.allOptionDefinitions) {
                val constant = property.toOptionConstant
                if (property.defaultValue !== null) {
                    members += property.toOptionDefault
                }
                if (property.lowerBound !== null) {
                    members += property.toOptionLowerBound
                }
                if (property.upperBound !== null) { 
                    members += property.toOptionUpperBound
                }
                members += constant
            }
            
            // 2. Private constants for required values of option dependencies
            for (property : bundle.members.allOptionDefinitions) {
                for (dependency : property.dependencies) {
                    if (dependency.value !== null) {
                        members += dependency.toDependencyValue
                    }
                }
            }
                        
            // 3. Implementation of ILayoutMetaDataProvider#apply(Registry)
            members += bundle.toMethod('apply', typeRef(void)) [
                parameters += bundle.toParameter('registry', typeRef("org.eclipse.elk.core.data.ILayoutMetaDataProvider.Registry"))
                body = '''
                    «registerLayoutOptions(bundle)»
                    «registerLayoutCategories(bundle)»
                    «registerLayoutAlgorithms(bundle)»
                '''
            ]
        ]
        
        // Infer code for the declared layout algorithms
        for (algorithm : bundle.members.filter(MdAlgorithm)) {
            inferAlgorithm(model, algorithm, acceptor, isPreIndexingPhase)
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Option Constants
    
    private def toOptionConstant(MdOption option) {
        return option.toField(
            option.constantName,
            typeRef(IProperty, option.type ?: typeRef(Object))
        ) [
            visibility = JvmVisibility.PUBLIC
            static = true
            final = true
            initializer = '''
                new «Property»<«option.type.asWrapperTypeIfPrimitive ?: typeRef(Object)»>(
                        «option.qualifiedName.toCodeString»«IF option.hasDefaultOrBounds»,
                        «IF option.defaultValue !== null»«option.defaultConstantName»«ELSE»null«ENDIF»,
                        «IF option.lowerBound !== null»«option.lowerBoundConstantName»«ELSE»null«ENDIF»,
                        «IF option.upperBound !== null»«option.upperBoundConstantName»«ELSE»null«ENDIF»«ENDIF»)'''
            if (option.deprecated) {
                annotations += annotationRef(Deprecated)
                deprecated = true
            }
            documentation = option.description.trimLines
        ]
    }
    
    private def hasDefaultOrBounds(MdOption option) {
        return option.defaultValue !== null || option.lowerBound !== null || option.upperBound !== null
    }
    
    private def toOptionDefault(MdOption option) {
        val optionType = option.type.cloneWithProxies ?: typeRef(Object)
        return option.toField(option.defaultConstantName, optionType) [
            visibility = JvmVisibility.PRIVATE
            static = true
            final = true
            initializer = option.defaultValue
            documentation = '''Default value for {@link #«option.constantName»}.'''
        ]
    }
    
    private def toOptionLowerBound(MdOption option) {
        val optionType = typeRef(Comparable, option.type.asWrapperTypeIfPrimitive?.wildcardSuper ?: wildcard)
        return option.toField(option.lowerBoundConstantName, optionType) [
            visibility = JvmVisibility.PRIVATE
            static = true
            final = true
            initializer = option.lowerBound
            documentation = '''Lower bound value for {@link #«option.constantName»}.'''
        ]
    }
    
    private def toOptionUpperBound(MdOption option) {
        val optionType = typeRef(Comparable, option.type.asWrapperTypeIfPrimitive?.wildcardSuper ?: wildcard)
        return option.toField(option.upperBoundConstantName, optionType) [
            visibility = JvmVisibility.PRIVATE
            static = true
            final = true
            initializer = option.upperBound
            documentation = '''Upper bound value for {@link #«option.constantName»}.'''
        ]
    }
    
    private def toDependencyValue(MdOptionDependency dependency) {
        val source = dependency.eContainer as MdOption
        val optionType = dependency.target.type.cloneWithProxies ?: typeRef(Object)
        return dependency.toField(dependency.dependencyConstantName, optionType) [
            visibility = JvmVisibility.PRIVATE
            static = true
            final = true
            initializer = dependency.value
            documentation = '''Required value for dependency between {@link #«source.constantName»} and {@link #«dependency.target.constantName»}.'''
        ]
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Registration Code
    
    private def StringConcatenationClient registerLayoutOptions(MdBundle bundle) '''
        «FOR option : bundle.members.getAllOptionDefinitions»
            registry.register(new «LayoutOptionData».Builder()
                .id(«option.qualifiedName.toCodeString»)
                .group(«option.groups.map[name].join('.').toCodeString»)
                .name(«(option.label ?: option.name).shrinkWhiteSpace.toCodeString»)
                .description(«option.description.shrinkWhiteSpace.toCodeString»)
                «IF option.defaultValue !== null»
                    .defaultValue(«option.defaultConstantName»)
                «ENDIF»
                «IF option.lowerBound !== null»
                    .lowerBound(«option.lowerBoundConstantName»)
                «ENDIF»
                «IF option.upperBound !== null»
                    .upperBound(«option.upperBoundConstantName»)
                «ENDIF»
                .type(«LayoutOptionData».Type.«option.optionType»)
                .optionClass(«option.optionTypeClass».class)
                «IF !option.targets.empty»
                    .targets(«EnumSet».of(«FOR t : option.targets SEPARATOR ', '»«LayoutOptionData».Target.«t.toString.toUpperCase»«ENDFOR»))
                «ENDIF»
                «IF option.programmatic || option.output || option.global»
                    .visibility(«LayoutOptionData».Visibility.HIDDEN)
                «ELSEIF option.advanced»
                    .visibility(«LayoutOptionData».Visibility.ADVANCED)
                «ELSE»
                    .visibility(«LayoutOptionData».Visibility.VISIBLE)
                «ENDIF»
                «IF !option.legacyIds.empty»
                    .legacyIds(«option.legacyIds.map['"' + it + '"'].join(', ')»)
                «ENDIF»
                .create()
            );
            «FOR dependency : option.dependencies»
                registry.addDependency(
                    «option.qualifiedName.toCodeString»,
                    «dependency.target.qualifiedName.toCodeString»,
                    «IF dependency.value === null»
                        null
                    «ELSE»
                        «dependency.dependencyConstantName»
                    «ENDIF»
                );
            «ENDFOR»
        «ENDFOR»
    '''
    
    private def LayoutOptionData() {
        typeRef('org.eclipse.elk.core.data.LayoutOptionData')
    }
    
    private def StringConcatenationClient registerLayoutCategories(MdBundle bundle) '''
        «FOR category : bundle.members.filter(MdCategory)»
            registry.register(new «LayoutCategoryData».Builder()
                .id(«category.qualifiedName.toCodeString»)
                .name(«(category.label ?: category.name).shrinkWhiteSpace.toCodeString»)
                .description(«category.description.shrinkWhiteSpace.toCodeString»)
                .create()
            );
        «ENDFOR»
    '''
    
    private def LayoutCategoryData() {
        typeRef('org.eclipse.elk.core.data.LayoutCategoryData')
    }
    
    private def StringConcatenationClient registerLayoutAlgorithms(MdBundle bundle) '''
        «FOR MdAlgorithm algorithm : bundle.members.filter(MdAlgorithm)»
            new «algorithm.qualifiedTargetClass»().apply(registry);
        «ENDFOR»
    '''
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ALGORITHM METADATA PROVIDER CLASS GENERATION

    /**
     * Generates code for the given algorithm's metadata provider class.
     * 
     * @param model
     *            the model the algorithm was declared in. Used for the file header.
     * @param algorithm
     *            the algorithm to create one or more {@link JvmDeclaredType declared types} from.
     * @param acceptor
     *            each created {@link JvmDeclaredType type} without a container should be passed to the acceptor
     *            in order to get attached to the current resource. The acceptor's
     *            {@link IJvmDeclaredTypeAcceptor#accept(org.eclipse.xtext.common.types.JvmDeclaredType) accept(..)}
     *            method takes the constructed empty type for the pre-indexing phase. This one is further
     *            initialized in the indexing phase using passed closure.
     * @param isPreIndexingPhase
     *            whether the method is called in a pre-indexing phase, i.e.
     *            when the global index is not yet fully updated. You must not
     *            rely on linking using the index if isPreIndexingPhase is
     *            <code>true</code>.
     */
    private def void inferAlgorithm(MdModel model, MdAlgorithm algorithm, IJvmDeclaredTypeAcceptor acceptor,
        boolean isPreIndexingPhase) {
            
        // If the algorithm doesn't have a name yet, don't bother
        if (algorithm.name === null) {
            return
        }
        
        // Create a class for the general option definitions to be used by clients as well as for the registration
        // method that instantiates the ILayoutMetaDataProviders
        acceptor.accept(algorithm.toClass(algorithm.qualifiedTargetClass)) [
            // Implement the ILayoutMetaDataProvider interface and... BE DOCUMENTED!
            superTypes += typeRef("org.eclipse.elk.core.data.ILayoutMetaDataProvider")
            fileHeader = model.documentation
            documentation = algorithm.documentation
            
            // 1. Public constants for supported layout options
            members += algorithm.toAlgorithmId
            
            for (support : algorithm.supportedOptions) {
                if (support.value !== null) {
                    members += support.toSupportedOptionDefault
                }
                members += support.toSupportedOptionConstant
            }
            
            // 2. Create a dedicated factory for the algorithm's layout provider.
            //    The generic 'AlgorithmFactory' isn't used on purpose: it's not supported by elk.js
            members += algorithm.toClass(algorithm.toFactoryClass) [
                superTypes += typeRef("org.eclipse.elk.core.util.IFactory<AbstractLayoutProvider>")
                static = true
                documentation = '''
                Layouter-specific algorithm factory.'''
                
                // Implementation of AbstractLayoutProvider#create()
                members += algorithm.toMethod('create', typeRef('org.eclipse.elk.core.AbstractLayoutProvider')) [
                    visibility = JvmVisibility.PUBLIC
                    body = '''
                        AbstractLayoutProvider provider = new «algorithm.provider»();
                        provider.initialize("«algorithm.parameter»");
                        return provider;
                    '''
                ]
                
                // Implementation of AbstractLayoutProvider#destroy(AbstractLayoutProvider)
                members += algorithm.toMethod('destroy', typeRef(void)) [
                    parameters += algorithm.toParameter('obj', typeRef('org.eclipse.elk.core.AbstractLayoutProvider'))
                    body = '''
                        obj.dispose();
                    '''
                ]
                
            ]
            
            // 3. Implementation of ILayoutMetaDataProvider#apply(Registry)
            members += algorithm.toMethod('apply', typeRef(void)) [
                parameters += algorithm.toParameter('registry', typeRef("org.eclipse.elk.core.data.ILayoutMetaDataProvider.Registry"))
                body = '''
                    «registerLayoutAlgorithm(algorithm)»
                '''
            ]
        ]
    }
    
    private def toFactoryClass(MdAlgorithm algorithm) {
        return algorithm.name.toFirstUpper + 'Factory'
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Algorithm ID Constant
    
    private def toAlgorithmId(MdAlgorithm algorithm) {
        return algorithm.toField("ALGORITHM_ID", typeRef(String)) [
            visibility = JvmVisibility.PUBLIC
            static = true
            final = true
            initializer = '''«algorithm.qualifiedName.toCodeString»'''
            documentation = '''The id of the «algorithm.label» algorithm.'''
        ]
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Property Constants
    
    private def toSupportedOptionDefault(MdOptionSupport support) {
        val algorithm = support.eContainer as MdAlgorithm
        val optionType = support.option.type.cloneWithProxies ?: typeRef(Object)
        return support.toField(support.option.defaultConstantName, optionType) [
            visibility = JvmVisibility.PRIVATE
            static = true
            final = true
            initializer = support.value
            documentation = '''Default value for {@link #«support.option.constantName»} with algorithm "«algorithm.label ?: algorithm.name»".'''
        ]
    }
    
    private def toSupportedOptionConstant(MdOptionSupport support) {
        return support.toField(
            support.option.constantName,
            typeRef(IProperty, support.option.type ?: typeRef(Object))
        ) [
            visibility = JvmVisibility.PUBLIC
            static = true
            final = true
            initializer = 
                if (support.value === null) {
                    '''«typeRef(support.option.bundle.qualifiedTargetClass)».«support.option.constantName»'''
                } else {
                    '''new «Property»<«support.option.type.asWrapperTypeIfPrimitive ?: typeRef(Object)»>(
                            «typeRef(support.option.bundle.qualifiedTargetClass)».«support.option.constantName»,
                            «support.option.defaultConstantName»)'''
                }
            if (support.option.deprecated) {
                annotations += annotationRef(Deprecated)
                deprecated = true
            }
            documentation = 
                support.option.description.trimLines 
                + if (!support.documentation.nullOrEmpty) "\n<h3>Algorithm Specific Details</h3>\n" + support.documentation.trimLines else ""
           
        ]
    }
    
    
    ////////////////////////////////////////////////////////////////////////////
    // Registration Code
    
    private def StringConcatenationClient registerLayoutAlgorithm(MdAlgorithm algorithm) '''
        registry.register(new «LayoutAlgorithmData».Builder()
            .id(«algorithm.qualifiedName.toCodeString»)
            .name(«(algorithm.label ?: algorithm.name).shrinkWhiteSpace.toCodeString»)
            .description(«algorithm.description.shrinkWhiteSpace.toCodeString»)
            .providerFactory(new «algorithm.toFactoryClass»())
            «IF algorithm.category !== null»
                .category(«algorithm.category.qualifiedName.toCodeString»)
            «ENDIF»
            «IF algorithm.bundle !== null»
                .melkBundleName(«algorithm.bundle.label.toCodeString»)
            «ENDIF»
            «IF algorithm.model?.definingBundleId !== null»
                .definingBundleId(«algorithm.model.definingBundleId.toCodeString»)
            «ENDIF»
            «IF algorithm.previewImage !== null»
                .imagePath(«algorithm.previewImage.toCodeString»)
            «ENDIF»
            «IF !algorithm.supportedFeatures.empty»
                .supportedFeatures(«EnumSet».of(«FOR f : algorithm.supportedFeatures SEPARATOR ', '»«typeRef('org.eclipse.elk.graph.properties.GraphFeature')».«f.toString.toUpperCase»«ENDFOR»))
            «ENDIF»
            «IF algorithm.validator !== null»
                .validatorClass(«algorithm.validator».class)
            «ENDIF»
            .create()
        );
        «FOR support : algorithm.supportedOptions»
            registry.addOptionSupport(
                «algorithm.qualifiedName.toCodeString»,
                «support.option.qualifiedName.toCodeString»,
                «IF support.value === null»
                    «support.option.constantName».getDefault()
                «ELSE»
                    «support.option.defaultConstantName»
                «ENDIF»
            );
        «ENDFOR»
    '''
    
    private def LayoutAlgorithmData() {
        typeRef('org.eclipse.elk.core.data.LayoutAlgorithmData')
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // UTILITY METHODS
    
    private def Iterable<MdOption> getAllOptionDefinitions(Iterable<? extends MdBundleMember> elements) {
        Iterables.concat(
            elements.filter(MdOption),
            elements.filter(MdGroup)
                    .map[it.children.getAllOptionDefinitions].flatten)
    }
    
    private def String getQualifiedTargetClass(MdBundle bundle) {
        val model = bundle.eContainer as MdModel
        val bundleClass = bundle.targetClass ?: 'Metadata'
        return model.name + '.' + bundleClass
    }
    
    private def String getQualifiedTargetClass(MdAlgorithm algorithm) {
        val model = algorithm.bundle.eContainer as MdModel
        val algorithmClass = algorithm.targetClass ?: algorithm.name.toFirstUpper + 'Metadata'
        return model.name + '.' + algorithmClass
    }
    
    private def Iterable<MdGroup> getGroups(MdBundleMember member) {
        val groups = new LinkedList
        var group = member.eContainer
        while (group instanceof MdGroup) {
            groups.addFirst(group)
            group = group.eContainer
        }
        groups
    }
    
    private def String getOptionType(MdOption option) {
        val jvmType = option.type?.type
        switch jvmType {
            
            JvmPrimitiveType: switch jvmType.identifier {
                case boolean.name:  return "BOOLEAN"
                case int.name:      return "INT"
                case double.name:    return "DOUBLE"
                // TODO it may be better to prevent float from the start
                case float.name:   return "DOUBLE"
            }
            
            JvmGenericType: switch jvmType.identifier {
                case Boolean.canonicalName:   return "BOOLEAN"
                case Integer.canonicalName:   return "INT"
                case Double.canonicalName:     return "DOUBLE"
                case Float.canonicalName:    return "DOUBLE"
                case String.canonicalName:    return "STRING"
                case EnumSet.canonicalName:   return "ENUMSET"
//                case jvmType.hasSupertype(IDataObject): return "OBJECT"
                default:                      return "OBJECT"
            }
            
            JvmEnumerationType: return "ENUM"
            
        }
        return "UNDEFINED";
    }  
    
    private def JvmTypeReference getOptionTypeClass(MdOption property) {
        if (property.type !== null) {
            // There are a few cases we need to catch
            if (property.type.type instanceof JvmPrimitiveType) {
                val primitiveType = property.type.type as JvmPrimitiveType;
                
                // Float and long need to be mapped to Double and Integer, respectively
                if (primitiveType.simpleName == "float") {
                    return typeRef(Double);
                } else if (primitiveType.simpleName == "long") {
                    return typeRef(Integer);
                }
            } else if (property.type.type instanceof JvmGenericType) {
                val genericType = property.type.type as JvmGenericType;
                
                // Double and Float... same as above. Also, with EnumSets we want to have the type of the set's
                // content as the option type, not the EnumSet type itself. Finally, with Lists, we need to drop the
                // generic type argument
                if (genericType.identifier == "java.lang.Float") {
                    return typeRef(Double);
                } else if (genericType.identifier == "java.lang.Long") {
                    return typeRef(Integer);
                } else if (genericType.identifier == "java.util.EnumSet") {
                    val outer = property.type as JvmParameterizedTypeReference
                    return outer.arguments.head.cloneWithProxies
                } else if (genericType.identifier == "java.util.List") {
                    return typeRef(List);
                }
            }
            
            // Standard case if we didn't run into the special cases
            property.type.cloneWithProxies.asWrapperTypeIfPrimitive
        } else {
            return typeRef(Void);
        }
    }
    
    private def MdBundle getBundle(MdBundleMember member) {
        var parent = member.eContainer
        while (!(parent instanceof MdBundle)) {
            parent = parent.eContainer
        }
        return parent as MdBundle
    }
    
    private def MdModel getModel(MdBundleMember member) {
        var parent = member.eContainer
        while (!(parent instanceof MdModel)) {
            parent = parent.eContainer
        }
        return parent as MdModel
    }
    
    private def String getDefiningBundleId(MdModel model) {
        val uri = model.eResource.URI
        val srcIndex = uri.segments.indexOf('src')
        if (srcIndex >= 1) {
            return uri.segments.get(srcIndex - 1)
        }
        return null
    }
    
    ////////////////////////////////////////////////////////////////////////////
    // Constant Names
    
    private def String getQualifiedName(MdBundleMember member) {
        val bundle = member.bundle
        val model = bundle.eContainer as MdModel
        var prefix = bundle.idPrefix ?: model.name
        // if the last part of the idPrefix would be repeated by the 
        // member's name, we skip it
        if (member.name == prefix.substring(prefix.lastIndexOf('.') + 1, prefix.length)) {
           prefix = prefix.substring(0, prefix.lastIndexOf('.')) 
        }
        return prefix
               + (if (member.groups.empty) '' else '.')
               + member.groups.map[it.name].join('.') 
               + '.' + member.name
    }
    
    private def String getConstantName(MdBundleMember member) {
        val name = member.name
        if (name !== null) {
            val result = new StringBuilder
            result.append(member.groups.map[it.name.toUpperCaseWithUnderscores].join('_'))
            if (result.length > 0) result.append('_')
            result.append(name.toUpperCaseWithUnderscores)
            return result.toString
        }
    }
    
    private def String toUpperCaseWithUnderscores(String str) {
        val result = new StringBuilder
        for (var i = 0; i < str.length; i++) {
            val c = str.charAt(i)
            if (Character.isUpperCase(c) && i > 0)
                result.append('_')
            result.append(Character.toUpperCase(c))
        }
        return result.toString
    }
    
    private def getDefaultConstantName(MdOption option) {
        option.constantName + '_DEFAULT'
    }
    
    private def getLowerBoundConstantName(MdOption option) {
        option.constantName + '_LOWER_BOUND'
    }
    
    private def getUpperBoundConstantName(MdOption option) {
        option.constantName + '_UPPER_BOUND'
    }
    
    private def getDependencyConstantName(MdOptionDependency dependency) {
        val option = dependency.eContainer as MdOption
        option.constantName + '_DEP_' + dependency.target.constantName + '_' + option.dependencies.indexOf(dependency)
    }
    
    private def String toCodeString(String s) {
        if (s === null)
            return 'null'
        else
            return '''"«Strings.convertToJavaString(s)»"'''
    }
    
    private def String shrinkWhiteSpace(String s) {
        if (s === null)
            return null
        val result = new StringBuilder
        var shrink = true
        for (var i = 0; i < s.length; i++) {
            val c = s.charAt(i)
            if (Character.isWhitespace(c)) {
                if (!shrink)
                    result.append(' ')
                shrink = true
            } else {
                result.append(c)
                shrink = false
            }
        }
        if (result.length > 0 && Character.isWhitespace(result.charAt(result.length - 1)))
            result.deleteCharAt(result.length - 1)
        return result.toString
    }
    
    private def String trimLines(String s) {
        if (s === null)
            null
        else
            s.split('\r?\n').map[trim].join('\n')
    }
    
}
