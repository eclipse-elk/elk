/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
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
package org.eclipse.elk.graph.text.validation

import com.google.inject.Inject
import java.util.EnumSet
import java.util.Map
import org.eclipse.elk.core.validation.LayoutOptionValidator
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.EMapPropertyHolder
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.properties.IPropertyValueProxy
import org.eclipse.elk.graph.util.ElkGraphUtil
import org.eclipse.emf.ecore.EObject
import org.eclipse.xtext.validation.Check

import static org.eclipse.elk.graph.ElkGraphPackage.Literals.*
import static org.eclipse.elk.graph.text.validation.IssueCodes.*

import static extension org.eclipse.elk.graph.text.ElkGraphTextUtil.*
import static extension org.eclipse.xtext.EcoreUtil2.*

/**
 * Custom validation rules for the ElkGraph language.
 */
class ElkGraphValidator extends AbstractElkGraphValidator {
    
    @Inject
    LayoutOptionValidator layoutOptionValidator
	
	@Check
	def void checkPropertyKey(ElkPropertyToValueMapEntryImpl entry) {
	   val optionData = LayoutMetaDataService.instance.getOptionData(entry.key.id)
	   if (optionData === null) {
	       warning("'" + entry.key.id + "' is not a layout option.", ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY)
	   }  
	}
	
	@Check
	def void checkPropertyValue(ElkPropertyToValueMapEntryImpl entry) {
	    val option = entry.key.toLayoutOption
	    if (option !== null) {
	        val container = entry.eContainer
	        if (container instanceof ElkGraphElement) {
	            switch container {
                    ElkNode:
                        checkOptionTarget(option, LayoutOptionData.Target.NODES, LayoutOptionData.Target.PARENTS)
                    ElkEdge:
                        checkOptionTarget(option, LayoutOptionData.Target.EDGES)
                    ElkPort:
                        checkOptionTarget(option, LayoutOptionData.Target.PORTS)
                    ElkLabel:
                        checkOptionTarget(option, LayoutOptionData.Target.LABELS)
                }
                if (container instanceof ElkNode && option.targets.contains(LayoutOptionData.Target.NODES))
                    checkAlgorithmSupport(option, (container as ElkNode).parent)
                else
                    checkAlgorithmSupport(option, container)
	        }
    	    var value = entry.value
    	    if (value instanceof IPropertyValueProxy) {
    	        value = value.resolveValue(option)
    	        if (value === null) {
        	        switch option.type {
        	            case STRING:
        	                expectPropertyType(String)
        	            case BOOLEAN:
                            expectPropertyType(Boolean)
        	            case INT:
                            expectPropertyType(Integer)
        	            case DOUBLE:
                            expectPropertyType(Double)
        	            case ENUMSET:
                            expectPropertyType(EnumSet)
        	            default:
        	                if (option.optionClass !== null)
        	                    expectPropertyType(option.optionClass)
        	        }
    	        } else {
    	            entry.value = value
    	        }
    	    }
    	    if (value !== null) {
                val issues = layoutOptionValidator.checkProperty(option, value,
                       entry.getContainerOfType(ElkGraphElement))
                for (issue : issues) {
                    switch issue.severity {
                        case ERROR:
                            error(issue.message, ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE)
                        case WARNING:
                            warning(issue.message, ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE)
                    }
                }
                if (CoreOptions.ALGORITHM == option) {
                    if (LayoutMetaDataService.instance.getAlgorithmDataBySuffix(value as String) === null)
                       error("No layout algorithm with identifier '" + value + "' can be found.",
                           ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE)
                }
            }
	    }
	}
    
    private def toLayoutOption(IProperty<?> property) {
        if (property instanceof LayoutOptionData)
            return property
        else if (property !== null)
            return LayoutMetaDataService.instance.getOptionData(property.id)
    }
    
    private def void checkOptionTarget(LayoutOptionData option, LayoutOptionData.Target... targetTypes) {
        if (!targetTypes.exists[option.targets.contains(it)])
            warning("The layout option '" + option.id + "' is not applicable to " + targetTypes.head.toString.toLowerCase + '.',
                ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY, OPTION_NOT_APPLICABLE)
    }
    
    private def void checkAlgorithmSupport(LayoutOptionData option, ElkGraphElement element) {
        if (option != CoreOptions.ALGORITHM) {
            var algorithm = element.algorithm
            if (algorithm !== null && !algorithm.knowsOption(option)) {
                // Maybe the option is handled by one of the parent nodes
                var foundMatch = false
                var parent = element.getContainerOfType(ElkNode)
                while (parent !== null && !foundMatch) {
                    algorithm = parent.algorithm
                    foundMatch = algorithm !== null && algorithm.knowsOption(option)
                    parent = parent.parent
                }
                if (!foundMatch)
                    warning("The algorithm '" + algorithm.id + "' does not support the option '" + option.id + "'.",
                        ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY)
            }
        }
    }
	
	private def void expectPropertyType(Class<?> type) {
       error("Expected value of type " + type.simpleName + ".", ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE)
	}
	
	@Check
	def void checkUniqueProperties(EMapPropertyHolder propertyHolder) {
	    val usedProperties = newHashMap
        for (entry : propertyHolder.properties) {
            val property = entry.key
            if (usedProperties.containsKey(property)) {
                propertyAlreadyAssigned(entry)
                val other = usedProperties.get(property)
                if (other !== null) {
                   propertyAlreadyAssigned(other)
                   usedProperties.put(property, null)
               }
            } else {
                usedProperties.put(property, entry)
            }
        }
	}
	
	private def void propertyAlreadyAssigned(Map.Entry<IProperty<?>, Object> entry) {
	    if (entry instanceof ElkPropertyToValueMapEntryImpl)
	        error("Property is already assigned.", entry, ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY)
	}
	
	@Check
	def void checkUniqueNames(ElkGraphElement element) {
	    val usedNames = newHashMap
	    for (object : element.eContents) {
	        val name = switch object {
	            ElkGraphElement: object.identifier
	            ElkEdgeSection: object.identifier
	        }
	        if (name !== null) {
	            if (usedNames.containsKey(name)) {
	                nameAlreadyUsed(object)
	                val other = usedNames.get(name)
	                if (other !== null) {
	                   nameAlreadyUsed(other)
	                   usedNames.put(name, null)
                   }
	            } else {
	                usedNames.put(name, object)
	            }
	        }
	    }
	}
	
	private def void nameAlreadyUsed(EObject object) {
	    val feature = switch object {
            ElkGraphElement: ELK_GRAPH_ELEMENT__IDENTIFIER
            ElkEdgeSection: ELK_EDGE_SECTION__IDENTIFIER
        }
	    error("Identifier is already used.", object, feature)
	}
	
	@Check
	def void checkEdgeContainer(ElkEdge edge) {
	    if (!(edge.sources.empty && edge.targets.empty)) {
	        val bestContainer = ElkGraphUtil.findBestEdgeContainment(edge)
	        if (bestContainer !== null && bestContainer != edge.containingNode) {
	            if (bestContainer.parent === null)
	                warning("This edge should be declared in the root node of this graph.", null)
	            else
	                warning("This edge should be declared in node " + bestContainer.identifier + ".", null)
	        }
        }
	}
	
	@Check
	def void checkEdgeSection(ElkEdgeSection edgeSection) {
	    val edge = edgeSection.parent
	    val incomingShape = edgeSection.incomingShape
        if (incomingShape !== null && !incomingShape.eIsProxy) {
            if (!edge.sources.contains(incomingShape))
	            error("The " + incomingShape.eClass.name + " " + incomingShape.identifier + " is not a source of this edge.",
	                ELK_EDGE_SECTION__INCOMING_SHAPE)
            if (!edgeSection.incomingSections.empty)
                error("An edge section cannot be connected to an " + incomingShape.eClass.name + " and other sections at the same time.",
                    ELK_EDGE_SECTION__INCOMING_SHAPE)
        }
        val outgoingShape = edgeSection.outgoingShape
        if (outgoingShape !== null && !outgoingShape.eIsProxy) {
            if (!edge.targets.contains(outgoingShape))
	            error("The " + outgoingShape.eClass.name + " " + outgoingShape.identifier + " is not a target of this edge.",
	                ELK_EDGE_SECTION__OUTGOING_SHAPE)
            if (!edgeSection.outgoingSections.empty)
                error("An edge section cannot be connected to an " + outgoingShape.eClass.name + " and other sections at the same time.",
                    ELK_EDGE_SECTION__OUTGOING_SHAPE)
        }
	}
	
}
