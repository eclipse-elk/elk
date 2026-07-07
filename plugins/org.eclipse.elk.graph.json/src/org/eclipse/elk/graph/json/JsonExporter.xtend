/**
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.json

import com.google.common.base.Splitter
import com.google.common.collect.BiMap
import com.google.common.collect.HashBiMap
import com.google.common.collect.Maps
import java.util.Collections
import java.util.Iterator
import java.util.Map
import java.util.Random
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.ElkShape
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.properties.IPropertyHolder

/**
 * Exporter from elk graph to json.
 */
final class JsonExporter {
    
    extension JsonAdapter = new JsonAdapter

    val BiMap<ElkNode, String> nodeIdMap = HashBiMap.create()
    val BiMap<ElkPort, String> portIdMap = HashBiMap.create()
    val BiMap<ElkEdge, String> edgeIdMap = HashBiMap.create()
    val BiMap<ElkEdgeSection, String> edgeSectionIdMap = HashBiMap.create()

    val Map<ElkNode, Object> nodeJsonMap = Maps.newHashMap
    val Map<ElkPort, Object> portJsonMap = Maps.newHashMap
    val Map<ElkEdge, Object> edgeJsonMap = Maps.newHashMap
    val Map<ElkEdgeSection, Object> edgeSectionJsonMap = Maps.newHashMap

    var nodeIdCounter = 0
    var portIdCounter = 0
    var edgeIdCounter = 0
    var edgeSectionIdCounter = 0
    
    // configuration
    var omitZeroPos = true
    var omitZeroDim = true
    var omitLayout = false
    var shortLayoutOptionKeys = true
    var omitUnknownLayoutOptions = true
    
    new () { }

    def setOptions(boolean omitZeroPos, boolean omitZeroDim, boolean omitLayout, 
        boolean shortLayoutOptionKeys, boolean omitUnknownLayoutOptions) {
        this.omitZeroPos = omitZeroPos
        this.omitZeroDim = omitZeroDim
        this.omitLayout = omitLayout
        this.shortLayoutOptionKeys = shortLayoutOptionKeys
        this.omitUnknownLayoutOptions = omitUnknownLayoutOptions
    }

    def export(ElkNode root) {
        init

        // create a tmp array
        val jsonArray = newJsonArray
        // transform the root node
        root.transformNode(jsonArray)

        // edges (important that all nodes are already transformed!)
        root.transformEdges

        // retrieve the result from the tmp array
        val jsonGraph = jsonArray.optJSONObject(0)
        return jsonGraph
    }

    private def init() {
        nodeIdMap.clear
        portIdMap.clear
        edgeIdMap.clear
        edgeSectionIdMap.clear
        nodeJsonMap.clear
        portJsonMap.clear
        edgeJsonMap.clear
        edgeSectionJsonMap.clear
        nodeIdCounter = 0
        portIdCounter = 0
        edgeIdCounter = 0
        edgeSectionIdCounter = 0
    }

    private def void transformNode(ElkNode node, Object arrayA) {
        // create the node and add it to the parent
        val jsonObj = node.createAndRegister
        val array = arrayA.toJsonArray
        array.addJsonArr(jsonObj)

        // labels
        if (!node.labels.nullOrEmpty) {
            val labels = newJsonArray
            jsonObj.addJsonObj("labels", labels)
            node.labels.forEach[ it.transformLabel(labels) ]
        }

        // ports
        if (!node.ports.nullOrEmpty) {
            val ports = newJsonArray
            jsonObj.addJsonObj("ports", ports)
            node.ports.forEach[ it.transformPort(ports) ]
        }

        // children
        if (!node.children.nullOrEmpty) {
            val children = newJsonArray
            jsonObj.addJsonObj("children", children)
            node.children.forEach[ it.transformNode(children) ]
        }
        
        // properties
        node.transformProperties(jsonObj)
        node.transformIndividualSpacings(jsonObj)
        node.transferShapeLayout(jsonObj)
    }

    private def void transformPort(ElkPort port, Object arrayA) {
        val jsonObj = port.createAndRegister
        val array = arrayA.toJsonArray
        array.addJsonArr(jsonObj)

        // labels
        if (!port.labels.nullOrEmpty) {
            val labels = newJsonArray
            jsonObj.addJsonObj("labels", labels)
            port.labels.forEach[ it.transformLabel(labels) ]
        }

        // properties and things
        port.transformProperties(jsonObj)
        port.transferShapeLayout(jsonObj)
    }

    private def void transformEdges(ElkNode node) {
        
        if (!node.containedEdges.nullOrEmpty) {
            val jsonObj = nodeJsonMap.get(node).toJsonObject
            val edges = newJsonArray
            jsonObj.addJsonObj("edges", edges)
            node.containedEdges.forEach[ it.transformEdge(edges) ]
        }

        // edges of children
        node.children.forEach[ it.transformEdges ]
    }

    private def void transformEdge(ElkEdge edge, Object arrayA) {
        val jsonObj = edge.createAndRegister
        val array = arrayA.toJsonArray
        array.addJsonArr(jsonObj)

        // connection points
        val sources = newJsonArray
        edge.sources.forEach[ s | 
            var source = portIdMap.get(s)
            if (source === null) {
                source = nodeIdMap.get(s)
            }
            if (source === null) {
                 formatError("Unknown edge source: " + s)
            }
            sources.addJsonArr(source.toJson)
        ]
        jsonObj.addJsonObj("sources", sources)
        
        val targets = newJsonArray
        edge.targets.forEach[ t | 
            var target = portIdMap.get(t)
            if (target === null) {
                target = nodeIdMap.get(t)
            }
            if (target === null) {
                formatError("Unknown edge target: " + target)
            }
            targets.addJsonArr(target.toJson)
        ]
        jsonObj.addJsonObj("targets", targets)

        // labels
        if (!edge.labels.nullOrEmpty) {
            val labels = newJsonArray
            jsonObj.addJsonObj("labels", labels)
            edge.labels.forEach[ it.transformLabel(labels) ]
        }
        
        // sections
        if (!omitLayout && !edge.sections.nullOrEmpty) {
            val sections = newJsonArray
            jsonObj.addJsonObj("sections", sections)
            edge.sections.forEach[ it.transformSection(sections) ]
        }

        // transfer junction points, if existent
        // make sure not to initialize an empty set of junction points by accident (#559)
        // when immediately using 'getProperty'
        if (!omitLayout && edge.hasProperty(CoreOptions.JUNCTION_POINTS)) {
            val jps = edge.getProperty(CoreOptions.JUNCTION_POINTS)
            if (!omitLayout && !jps.nullOrEmpty) {
                val jsonJPs = newJsonArray
                jps.forEach[ jp |
                    val jsonPnt = newJsonObject
                    jsonPnt.addJsonObj("x", jp.x)
                    jsonPnt.addJsonObj("y", jp.y)
                    jsonJPs.addJsonArr(jsonPnt)
                ]
                jsonObj.addJsonObj("junctionPoints", jsonJPs)
            }
        }

        // properties        
        edge.transformProperties(jsonObj)
    }
    
    private def void transformSection(ElkEdgeSection section, Object sectionsA) {
        val jsonObj = section.createAndRegister
        val sections = sectionsA.toJsonArray
        sections.addJsonArr(jsonObj)
        
        // Start Point
        val startPoint = newJsonObject
        startPoint.addProperty("x", section.startX)
        startPoint.addProperty("y", section.startY)
        jsonObj.addJsonObj("startPoint", startPoint)
        
        // End Point
        val endPoint = newJsonObject
        endPoint.addProperty("x", section.endX)
        endPoint.addProperty("y", section.endY)
        jsonObj.addJsonObj("endPoint", endPoint)
        
        // Bend Points
        if (!omitLayout && !section.bendPoints.nullOrEmpty) {
            val bendPoints = newJsonArray
                section.bendPoints.forEach [ pnt |
                    val jsonPnt = newJsonObject
                    jsonPnt.addProperty("x", pnt.x)
                    jsonPnt.addProperty("y", pnt.y)
                    bendPoints.addJsonArr(jsonPnt)
                ]
            jsonObj.addJsonObj("bendPoints", bendPoints)
        }
        
        // Incoming shape
        if (section.incomingShape !== null) {
            jsonObj.addProperty("incomingShape", idByElement(section.incomingShape))
        }
        
        // Outgoing shape
        if (section.outgoingShape !== null) {
            jsonObj.addProperty("outgoingShape", idByElement(section.outgoingShape))
        }
        
        // Incoming sections
        if (!section.incomingSections.nullOrEmpty) {
            val incomingSections = newJsonArray
            section.incomingSections.forEach [ sec | incomingSections.addJsonArr(idByElement(sec).toJson) ]
            jsonObj.addJsonObj("incomingSections", incomingSections)
        }
        
        // Outgoing sections
        if (!section.outgoingSections.nullOrEmpty) {
            val outgoingSections = newJsonArray
            section.outgoingSections.forEach [ sec | outgoingSections.addJsonArr(idByElement(sec).toJson) ]
            jsonObj.addJsonObj("outgoingSections", outgoingSections)
        }
        
        section.transformProperties(jsonObj);
    }

    private def void transformLabel(ElkLabel label, Object array) {
        val jsonLabel = newJsonObject
        jsonLabel.addProperty("text", label.text)
        if (!label.identifier.nullOrEmpty) {
            jsonLabel.addProperty("id", label.identifier)
        }
        array.toJsonArray.addJsonArr(jsonLabel)

        // properties and things
        label.transformProperties(jsonLabel)
        label.transferShapeLayout(jsonLabel)
    }

    private def void transformProperties(IPropertyHolder holder, Object parentA) {
        // skip if empty
        if (holder === null || holder.allProperties === null || holder.allProperties.empty) {
            return
        }
        
        val jsonProps = newJsonObject
        val parent = parentA.toJsonObject
        parent.addJsonObj("layoutOptions", jsonProps)
        holder.getAllProperties.entrySet
            .filter[ key !== null ]
            .filter[ key != CoreOptions.SPACING_INDIVIDUAL ]
            .forEach [ p |
                if (!omitUnknownLayoutOptions || p.key.isKnown) {
                    var key = if (shortLayoutOptionKeys) p.key.id.shortOptionKey else p.key.id
                    jsonProps.addProperty(key, p.value.toString)                
                }
            ]
    }
    
    private def void transformIndividualSpacings(IPropertyHolder holder, Object parentA) {
        // skip if empty
        if (holder === null || !holder.hasProperty(CoreOptions.SPACING_INDIVIDUAL)) {
            return
        }
        val individualSpacings = holder.getProperty(CoreOptions.SPACING_INDIVIDUAL)
        if (individualSpacings.allProperties === null || individualSpacings.allProperties.empty) {
            return;
        }
        
        val jsonProps = newJsonObject
        val parent = parentA.toJsonObject
        parent.addJsonObj("individualSpacings", jsonProps)
        individualSpacings.allProperties.entrySet.filter[it.key !== null].forEach [ p |
            if (!omitUnknownLayoutOptions || p.key.isKnown) {
                var key = if (shortLayoutOptionKeys) p.key.id.shortOptionKey else p.key.id
                jsonProps.addProperty(key, p.value.toString)                
            }
        ]
    }

    private def transferShapeLayout(ElkShape shape, Object jsonObjA) {

        val jsonObj = jsonObjA.toJsonObject
        // position
        if (!omitLayout) {
            if (!omitZeroPos) {
                jsonObj.addProperty("x", shape.x)
                jsonObj.addProperty("y", shape.y)  
            } else {
                // non-equality with double is fine here
                if (shape.x != 0.0) jsonObj.addProperty("x", shape.x)
                if (shape.y != 0.0) jsonObj.addProperty("y", shape.y)
            }
        }
        // dimension
        if (!omitZeroDim) {
            jsonObj.addProperty("width", shape.width)
            jsonObj.addProperty("height", shape.height)    
        } else {
            if (shape.width != 0.0) jsonObj.addProperty("width", shape.width)
            if (shape.height != 0.0) jsonObj.addProperty("height", shape.height)
        }
    }

    /* ---------------------------------------------------------------------------
     *   Convenience methods
     */
    private def createAndRegister(ElkNode node) {
        val obj = newJsonObject
        var id = node.identifier 
        if (id === null) {
            id = "n" + nodeIdCounter
            nodeIdCounter = nodeIdCounter + 1              
        }
        id = id.assertUnique(nodeIdMap.inverse)
        
        obj.addProperty("id", id)
        nodeIdMap.put(node, id)
        nodeJsonMap.put(node, obj)

        return obj
    }

    private def createAndRegister(ElkPort port) {
        val obj = newJsonObject
        var id = port.identifier 
        if (id === null) {
            id = "p" + portIdCounter
            portIdCounter = portIdCounter + 1    
        }
        id = id.assertUnique(portIdMap.inverse)

        obj.addProperty("id", id)
        portIdMap.put(port, id)
        portJsonMap.put(port, obj)

        return obj
    }

    private def createAndRegister(ElkEdge edge) {
        val obj = newJsonObject
        var id = edge.identifier
        if (id === null) {
            id = "e" + edgeIdCounter
            edgeIdCounter = edgeIdCounter + 1
        }
        id = id.assertUnique(edgeIdMap.inverse)

        obj.addProperty("id", id)
        edgeIdMap.put(edge, id)
        edgeJsonMap.put(edge, obj)

        return obj
    }

    private def createAndRegister(ElkEdgeSection section) {
        val obj = newJsonObject
        var id = section.identifier
        if (id === null) {
            id = "s" + edgeSectionIdCounter
            edgeSectionIdCounter = edgeSectionIdCounter + 1
        }
        id = id.assertUnique(edgeSectionIdMap.inverse)

        obj.addProperty("id", id)
        edgeSectionIdMap.put(section, id)
        edgeSectionJsonMap.put(section, obj)

        return obj
    }
    
    private def dispatch idByElement(ElkNode node) {
        return nodeIdMap.get(node);
    }
    
    private def dispatch idByElement(ElkPort port) {
        return portIdMap.get(port);
    }
    
    private def dispatch idByElement(ElkEdgeSection section) {
        return edgeSectionIdMap.get(section);
    }
    
    private def getShortOptionKey(String fullId) {
        val option = LayoutMetaDataService.instance.getOptionDataBySuffix(fullId)
        if (option === null) {
            // if the option is unknown, return the full id
            return fullId
        }
        val idSplit = Splitter.on('.').split(option.id)
        var foundMatch = false
        var i = idSplit.size - 1
        if (i >= 1 && option.group == idSplit.get(i - 1)) {
            i--
        }
        while (i >= 0 && !foundMatch) {
            val suffix = idSplit.drop(i)
            if (LayoutMetaDataService.instance.getOptionDataBySuffix(suffix.join('.')) !== null) {
                foundMatch = true
            } else {
                i--                
            }
        }
        if (foundMatch) {
            return idSplit.drop(i).join('.')
        } else {
            return option.id
        }
    }
        
    def <T> Iterator<T> emptyIfNull(Iterator<T> iterator) {
        if (iterator === null) {
            return Collections.emptyIterator
        } else {
            iterator
        }
    } 
    
    private def isKnown(IProperty<?> property) {
        return LayoutMetaDataService.instance.getOptionDataBySuffix(property.id) !== null
    }
    
    val RANDOM = new Random()
    private def String sixDigitRandomNumber() {
        return RANDOM.nextInt(1000000) + ""
    }
    
    private def String padZeroes(String s, int length) {
        var tmp = s
        while (tmp.length < length) {
            tmp = "0" + tmp
        }
        return tmp
    }
    
    private def String assertUnique(String id, Map<String, ?> map) {
        var tmp = id
        while (map.containsKey(tmp)) {
            tmp = id + "_g" + sixDigitRandomNumber().padZeroes(6)
        }
        return tmp
    }
}
