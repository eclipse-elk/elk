/**
 * Copyright (c) 2017, 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 */
package org.eclipse.elk.graph.json

import com.google.common.collect.BiMap
import com.google.common.collect.HashBiMap
import com.google.common.collect.HashMultimap
import com.google.common.collect.Maps
import com.google.common.collect.Multimap
import java.util.Map
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.core.options.EdgeCoords
import org.eclipse.elk.core.options.ShapeCoords
import org.eclipse.elk.core.util.IndividualSpacings
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.elk.graph.ElkShape
import org.eclipse.elk.graph.properties.IPropertyHolder
import org.eclipse.elk.graph.util.ElkGraphUtil

/**
 * Importer from json to elk graph. Internally it maintains a mapping that can be used to 
 * transfer any computed layout information to the original json. 
 * 
 * <h3>Implementation Hints</h3>
 * The implementation of the importer is kept free of any explicit json library. This is possible thanks 
 * to xtend's extensions methods. Wherever possible type-inference is used, if not the Object type is used.
 * The library-dependent code can be found in the {@link JsonAdapter} class. 
 * 
 * To get an id that must be specified and to preserve the id's type (string or int) use 'getId' 
 * (throws an exeption otherwise). If the result may be null and may always be a string
 * (e.g. when assembling the text for an exception) use 'getIdSave'. 
 */
final class JsonImporter {

    extension JsonAdapter = new JsonAdapter

    /* Id -> ElkGraph element maps
     * Id can be string or integer, thus {@link Object} is used. */
    val BiMap<Object, ElkNode> nodeIdMap = HashBiMap.create()
    val BiMap<Object, ElkPort> portIdMap = HashBiMap.create()
    val Map<Object, ElkEdge> edgeIdMap = Maps.newHashMap
    val BiMap<Object, ElkEdgeSection> edgeSectionIdMap = HashBiMap.create()

    /* ElkGraph element -> Json element maps */
    val BiMap<ElkNode, Object> nodeJsonMap = HashBiMap.create()
    val Map<ElkPort, Object> portJsonMap = Maps.newHashMap
    val Map<ElkEdge, Object> edgeJsonMap = Maps.newHashMap
    val Map<ElkEdgeSection, Object> edgeSectionJsonMap = Maps.newHashMap
    val Map<ElkLabel, Object> labelJsonMap = Maps.newHashMap
    
    /* Maps to help in adjusting coordinates */
    val Map<ElkEdge, ElkNode> edgeOriginalParentMap = Maps.newHashMap
    val Map<ElkShape, Double> globalXMap = Maps.newHashMap
    val Map<ElkShape, Double> globalYMap = Maps.newHashMap
    val Map<ElkGraphElement, ShapeCoords> shapeCoordsMap = Maps.newHashMap
    val Map<ElkGraphElement, EdgeCoords> edgeCoordsMap = Maps.newHashMap

    var Object inputModel

    /* ---------------------------------------------------------------------------
     *   JSON --> ElkGraph
     */
    /**
      * Main entry point for the json to ELK graph transformation. Runs through all elements
      * of the graph (nodes, ports, edges, edge sections) and creates correlating ELK graph elements.
      */
    def ElkNode transform(Object graph) {
        inputModel = graph
        clearMaps

        // transform the root node along with its children
        val root = graph.transformNode(null)
       
        // transform all edges
        graph.transformEdges

        // return the transformed ELK graph
        return root
    }

    def getInputModel() {
        return inputModel
    }

    private def clearMaps() {
        nodeIdMap.clear
        portIdMap.clear 
        edgeIdMap.clear
        edgeSectionIdMap.clear
        
        nodeJsonMap.clear
        portJsonMap.clear
        edgeJsonMap.clear
        edgeSectionJsonMap.clear
        labelJsonMap.clear
        
        edgeOriginalParentMap.clear
        globalXMap.clear
        globalYMap.clear
        shapeCoordsMap.clear
        edgeCoordsMap.clear
    }

    private def transformChildNodes(Object jsonNodeA, ElkNode parent) {
        val jsonNode = jsonNodeA.toJsonObject
        jsonNode.optJSONArray("children") => [ children |
            if (children !== null) {
                for (i : 0 ..< children.sizeJsonArr) {
                    children.optJSONObject(i)?.transformNode(parent)
                }
            }
        ]
    }

    private def ElkNode transformNode(Object jsonNode, ElkNode parent) {
        // create an ElkNode and add it to the parent
        val node = ElkGraphUtil.createNode(parent).register(jsonNode)
        node.identifier = jsonNode.toJsonObject.idSave 
        
        jsonNode.transformProperties(node)
        jsonNode.transformIndividualSpacings(node)
        jsonNode.transformShapeLayout(node)
        jsonNode.transformPorts(node)
        jsonNode.transformLabels(node)
        jsonNode.transformChildNodes(node)

        return node
    }

    private def void transformEdges(Object jsonObjA) {
        val jsonObj = jsonObjA.toJsonObject
        // the json object represents a node
        val node = nodeJsonMap.inverse.get(jsonObj);
        if (node === null) {
            throw formatError("Unable to find elk node for json object '" + jsonObj.idSave + "' Panic!")
        }
        
        // transform edges of the current hierarchy level
        jsonObj.optJSONArray("edges") => [ edges |
            if (edges !== null) {
                for (i : 0 ..< edges.sizeJsonArr) {
                    val edge = edges.optJSONObject(i)
                    if (edge.hasJsonObj("sources") || edge.hasJsonObj("targets")) {
                        edge.transformEdge(node)
                    } else {
                        edge.transformPrimitiveEdge(node)
                    }
                    ElkGraphUtil.updateContainment(edgeIdMap.get(edge.id))
                }
            }
        ]

        // transform the edges of all child nodes
        jsonObj.optJSONArray("children") => [ children |
            if (children !== null) {
                for (i : 0 ..< children.sizeJsonArr) {
                    children.optJSONObject(i)?.transformEdges
                }
            }
        ]
    }
    
    private def transformPrimitiveEdge(Object jsonObjA, ElkNode parent) {
        val jsonObj = jsonObjA.toJsonObject
        // Create ElkEdge
        val edge = ElkGraphUtil.createEdge(parent).register(jsonObj)
        edge.identifier = jsonObj.idSave

        // source
        val srcNode = nodeIdMap.get(jsonObj.getJsonObj("source").asId)
        val srcPort = portIdMap.get(jsonObj.getJsonObj("sourcePort")?.asId) // may be null

        if (srcNode === null) {
            throw formatError("An edge must have a source node (edge id: '" + jsonObj.id +"').") 
        }
        if (srcPort !== null && srcPort.parent != srcNode) {
            throw formatError("The source port of an edge must be a port of the edge's source node (edge id: '" 
                + jsonObj.idSave + "').")
        }
        
        edge.sources += srcPort ?: srcNode

        // target
        val tgtNode = nodeIdMap.get(jsonObj.getJsonObj("target").asId)
        val tgtPort = portIdMap.get(jsonObj.getJsonObj("targetPort")?.asId) // may be null
        
        if (tgtNode === null) {
            throw formatError("An edge must have a target node (edge id: '" + jsonObj.id +"').") 
        }
        
        if (tgtPort !== null && tgtPort.parent != tgtNode) {
            throw formatError("The target port of an edge must be a port of the edge's target node (edge id: '" 
                + jsonObj.idSave + "').")
        }
        
        edge.targets += tgtPort ?: tgtNode
        
        // check if ok
        if (edge.sources.empty || edge.targets.empty) {
            throw formatError("An edge must have at least one source and one target (edge id: '" 
                + jsonObj.idSave + "').")
        }
        
        jsonObj.transformProperties(edge)
        jsonObj.transformPrimitiveEdgeLayout(edge)
        jsonObj.transformLabels(edge)
    }

    private def transformPrimitiveEdgeLayout(Object jsonObjA, ElkEdge edge) {
        val jsonObj = jsonObjA.toJsonObject
        
        // if there's any edge route information available, turn it into an ElkSection
        if (jsonObj.hasJsonObj("sourcePoint") 
            || jsonObj.hasJsonObj("targetPoint") 
            || jsonObj.hasJsonObj("bendPoints")) { 

            val section = ElkGraphUtil.createEdgeSection(edge)
            // src
            jsonObj.optJSONObject("sourcePoint") => [ srcPnt |
                if (srcPnt !== null) {
                    srcPnt.optDouble("x") => [section.startX = it]
                    srcPnt.optDouble("y") => [section.startY = it]
                }
            ]
    
            // tgt
            jsonObj.optJSONObject("targetPoint") => [ tgtPnt |
                if (tgtPnt !== null) {
                    tgtPnt.optDouble("x") => [section.endX = it]
                    tgtPnt.optDouble("y") => [section.endY = it]
                }
            ]
    
            // bend points
            jsonObj.optJSONArray("bendPoints") => [ bendPoints |
                if (bendPoints !== null) {
                    for (i : 0 ..< bendPoints.sizeJsonArr) {
                        bendPoints.optJSONObject(i) => [ bendPoint |
                            ElkGraphUtil.createBendPoint(section, bendPoint.optDouble("x"), bendPoint.optDouble("y"));
                        ]
                    }
                }
            ]
        }
    }
    

    private def transformEdge(Object jsonObjA, ElkNode parent) {
        val jsonObj = jsonObjA.toJsonObject
        // Create ElkEdge
        val edge = ElkGraphUtil.createEdge(parent).register(jsonObj)
        edge.identifier = jsonObj.idSave

        // sources
        jsonObj.optJSONArray("sources") => [ sources |
            if (sources !== null) {
                for (i : 0 ..< sources.sizeJsonArr) {
                    val sourceElement = shapeById(sources.getJsonArr(i).asId);
                    if (sourceElement !== null) {
                        edge.sources += sourceElement;
                    }
                }
            }
        ]

        // targets
        jsonObj.optJSONArray("targets") => [ targets |
            if (targets !== null) {
                for (i : 0 ..< targets.sizeJsonArr) {
                    val targetElement = shapeById(targets.getJsonArr(i).asId);
                    if (targetElement !== null) {
                        edge.targets += targetElement;
                    }
                }
            }
        ]
        
        // check if ok
        if(edge.sources.empty || edge.targets.empty) {
            throw formatError("An edge must have at least one source and one target (edge id: '" 
                + jsonObj.idSave + "').")
        }
        
        // transform things
        jsonObj.transformProperties(edge)
        jsonObj.transformEdgeSections(edge)
        jsonObj.transformLabels(edge)
    }

    private def transformEdgeSections(Object jsonObjA, ElkEdge edge) {
        val jsonObj = jsonObjA.toJsonObject
        // While iterating over the edge's edge sections, we remember identifiers of the section's incoming and
        // outgoing edge sections. Those references, along with one special case for incoming and outgoing shapes,
        // are resolved later, after all sections have been transformed
        val Multimap<ElkEdgeSection, Object> incomingSectionIdentifiers = HashMultimap.create();
        val Multimap<ElkEdgeSection, Object> outgoingSectionIdentifiers = HashMultimap.create();
        
        jsonObj.optJSONArray("sections") => [ sections |
            if (sections !== null) {
                for (i : 0 ..< sections.sizeJsonArr) {
                    sections.optJSONObject(i) => [ jsonSection |
                        val elkSection = ElkGraphUtil.createEdgeSection(edge).register(jsonSection)
                        elkSection.identifier = jsonSection.idSave
                        
                        fillEdgeSectionCoordinates(jsonSection, elkSection);
                        
                        // Incoming and Outgoing shapes
                        jsonSection.optString("incomingShape") => [ jsonShapeId |
                            if (jsonShapeId !== null) {
                                elkSection.incomingShape = shapeById(jsonShapeId);
                            }
                        ]
                        
                        jsonSection.optString("outgoingShape") => [ jsonShapeId |
                            if (jsonShapeId !== null) {
                                elkSection.outgoingShape = shapeById(jsonShapeId);
                            }
                        ]
                        
                        // References to incoming and outgoing sections
                        jsonSection.optJSONArray("incomingSections") => [ jsonSectionIds |
                            if (jsonSectionIds !== null) {
                                for (j : 0 ..< jsonSectionIds.sizeJsonArr) {
                                    incomingSectionIdentifiers.put(elkSection, jsonSectionIds.getJsonArr(j).asId)
                                }
                            }
                        ]
                        
                        jsonSection.optJSONArray("outgoingSections") => [ jsonSectionIds |
                            if (jsonSectionIds !== null) {
                                for (j : 0 ..< jsonSectionIds.sizeJsonArr) {
                                    outgoingSectionIdentifiers.put(elkSection, jsonSectionIds.getJsonArr(j).asId)
                                }
                            }
                        ]
                    ]
                }
            }
        ]
        
        // Fill in references to incoming and outgoing sections
        for (section : incomingSectionIdentifiers.keySet) {
            for (id : incomingSectionIdentifiers.get(section)) {
                val referencedSection = edgeSectionIdMap.get(id);
                if (referencedSection !== null) {
                    section.incomingSections += referencedSection;
                } else {
                    throw formatError("Referenced edge section does not exist: " + id 
                        + " (edge id: '" + jsonObj.idSave + "').")
                }
            }
        }
        
        for (section : outgoingSectionIdentifiers.keySet) {
            for (id : outgoingSectionIdentifiers.get(section)) {
                val referencedSection = edgeSectionIdMap.get(id);
                if (referencedSection !== null) {
                    section.outgoingSections += referencedSection;
                } else {
                    throw formatError("Referenced edge section does not exist: " + id 
                        + " (edge id: '" + jsonObj.idSave + "').")
                }
            }
        }
        
        // Special case: if the edge has only a single source, a single target, and a single edge section which has
        // no incoming and outgoing shapes, set the incoming and outgoing shape to the source and target of the edge,
        // respectively
        if (edge.isConnected && !edge.isHyperedge && edge.sections.size == 1) {
            val section = edge.sections.get(0);
            if (section.incomingShape === null && section.outgoingShape === null) {
                section.incomingShape = edge.sources.get(0);
                section.outgoingShape = edge.targets.get(0);
            }
        }
    }
    
    private def fillEdgeSectionCoordinates(Object jsonObjA, ElkEdgeSection section) {
        val jsonObj = jsonObjA.toJsonObject
        jsonObj.optJSONObject("startPoint") => [ startPoint |
            if (startPoint !== null) {
                startPoint.optDouble("x") => [section.startX = it]
                startPoint.optDouble("y") => [section.startY = it]
            } else {
                throw formatError("All edge sections need a start point.")
            }
        ]
        
        jsonObj.optJSONObject("endPoint") => [ endPoint |
            if (endPoint !== null) {
                endPoint.optDouble("x") => [section.endX = it]
                endPoint.optDouble("y") => [section.endY = it]
            } else {
                throw formatError("All edge sections need an end point.")
            }
        ]
        
        jsonObj.optJSONArray("bendPoints") => [ bendPoints |
            if (bendPoints !== null) {
                for (i : 0 ..< bendPoints.sizeJsonArr) {
                    bendPoints.optJSONObject(i) => [ bendPoint |
                        ElkGraphUtil.createBendPoint(section, bendPoint.optDouble("x"), bendPoint.optDouble("y"));
                    ]
                }
            }
        ]
    }

    private def transformProperties(Object jsonObjA, IPropertyHolder layoutData) {
        val jsonObj = jsonObjA.toJsonObject
        
        var layoutOptions = jsonObj.optJSONObject("layoutOptions")
        if (layoutOptions === null) {
            // legacy support
            layoutOptions = jsonObj.optJSONObject("properties")
        }
        
        if (layoutOptions !== null) {
            val opts = layoutOptions
            opts?.keysJsonObj?.forEach[ k |
                val value = opts.getJsonObj(k)?.stringVal
                layoutData.setOption(k, value)
            ]
        }
    }
    
    private def transformIndividualSpacings(Object jsonObjA, IPropertyHolder layoutData) {
        val jsonObj = jsonObjA.toJsonObject
        
        val jsonIndividualSpacings = jsonObj.optJSONObject("individualSpacings")
        if (jsonIndividualSpacings !== null) {
            if (!layoutData.hasProperty(CoreOptions.SPACING_INDIVIDUAL)) {
                layoutData.setProperty(CoreOptions.SPACING_INDIVIDUAL, new IndividualSpacings())   
            }
            val individualSpacings = layoutData.getProperty(CoreOptions.SPACING_INDIVIDUAL);
            val opts = jsonIndividualSpacings
            opts?.keysJsonObj?.forEach[ k |
                val value = opts.getJsonObj(k)?.stringVal
                individualSpacings.setOption(k, value)
            ]
        }
    }
       
    private def setOption(IPropertyHolder e, String id, String value) {
        val optionData = LayoutMetaDataService.instance.getOptionDataBySuffix(id)
        if (optionData !== null) {
            val parsed = optionData.parseValue(value)
            if (parsed !== null) {
                e.setProperty(optionData, parsed)
            }
        }
    } 

    private def transformLabels(Object jsonObjA, ElkGraphElement element) {
        val jsonObj = jsonObjA.toJsonObject
        jsonObj.optJSONArray("labels") => [ labels |
            if (labels !== null) {
                for (i : 0 ..< labels.sizeJsonArr) {
                    val jsonLabel = labels.optJSONObject(i)
                    if (jsonLabel !== null) {
                        val label = ElkGraphUtil.createLabel(jsonLabel.optString("text"), element)
                        labelJsonMap.put(label, jsonLabel) 
                        if (jsonLabel.hasJsonObj("id")) {
                            label.identifier = jsonLabel.idSave
                        }
                        
                        jsonLabel.transformProperties(label)
                        jsonLabel.transformShapeLayout(label)
                    }
                }
            }
        ]
    }

    private def transformPorts(Object jsonObjA, ElkNode parent) {
        val jsonObj = jsonObjA.toJsonObject
        jsonObj.optJSONArray("ports") => [ ports |
            if (ports !== null) {
                for (i : 0 ..< ports.sizeJsonArr) {
                    ports.optJSONObject(i)?.transformPort(parent)
                }
            }
        ]
    }

    private def transformPort(Object jsonPort, ElkNode parent) {
        // create ElkPort
        val port = ElkGraphUtil.createPort(parent).register(jsonPort)
        port.identifier = jsonPort.toJsonObject.idSave

        // transform things
        jsonPort.transformProperties(port)
        jsonPort.transformShapeLayout(port)
        jsonPort.transformLabels(port)
    }
    
    private def transformShapeLayout(Object jsonObjA, ElkShape shape) {
        val jsonObj = jsonObjA.toJsonObject
        jsonObj.optDouble("x") => [shape.x = it.doubleValueValid]
        jsonObj.optDouble("y") => [shape.y = it.doubleValueValid]
        jsonObj.optDouble("width") => [shape.width = it.doubleValueValid]
        jsonObj.optDouble("height") => [shape.height = it.doubleValueValid]
    }
    
    private def double doubleValueValid(Double d) {
        if (d === null || d.infinite || d.naN) {
            return 0.0
        } else {
            return d.doubleValue
        }
    }
    
    private def shapeById(Object id) {
        val node = nodeIdMap.get(id)
        if (node !== null) {
            return node
        } 
        val port = portIdMap.get(id)
        if (port !== null) {
            return port
        }
        throw formatError("Referenced shape does not exist: " + id)
    }

    /* ---------------------------------------------------------------------------
     *   ElkGraph positions -> Json
     */
    /**
      * Transfer the layout back to the formerly imported graph, using {@link #transform(Object)}.
      */
    def transferLayout(ElkNode graph) {
        // First pass handles nodes and ports.
        ElkGraphUtil.propertiesSkippingIteratorFor(graph, true).forEach [ element |
            element.transferLayoutInt1
        ]
        
        // Second pass handles edges and labels.
        ElkGraphUtil.propertiesSkippingIteratorFor(graph, true).forEach [ element |
            element.transferLayoutInt2
        ]
    }

    private def dispatch transferLayoutInt1(ElkNode node) {
        val jsonObj = nodeJsonMap.get(node)
        if (jsonObj === null) {
            throw formatError("Node did not exist in input.")
        }
        node.recordGlobalCoords
        node.recordCoordinateModes
        node.transferShapeLayout(jsonObj)
    }

    private def dispatch transferLayoutInt1(ElkPort port) {
        val jsonObj = portJsonMap.get(port)
        if (jsonObj === null) {
            throw formatError("Port did not exist in input.")
        }
        port.recordGlobalCoords
        port.recordCoordinateModes
        port.transferShapeLayout(jsonObj)
    }

    private def dispatch transferLayoutInt2(ElkEdge edge) {
        val jsonObj = edgeJsonMap.get(edge).toJsonObject
        if (jsonObj === null) {
            throw formatError("Edge did not exist in input.")
        }
        
        edge.recordCoordinateModes
        
        val edgeId = jsonObj.id
                
        // what we need to transfer are the edge sections
        if (!edge.sections.nullOrEmpty) {
            val sections = newJsonArray
            edge.sections.forEach [ elkSection, i |
                
                // check if there's a corresponding json section
                var maybeSection = edgeSectionJsonMap.get(elkSection)
                if (maybeSection === null) {
                    // otherwise create a new one
                    maybeSection = newJsonObject 
                    // Id, just enumerate the sections per edge
                    maybeSection.toJsonObject.addJsonObj("id", edgeId + "_s" + i)
                }
                val jsonSection = maybeSection.toJsonObject
                sections.addJsonArr(jsonSection)
                
                // Start Point
                val startPoint = newJsonObject
                startPoint.addJsonObj("x", edge.adjustX(elkSection.startX))
                startPoint.addJsonObj("y", edge.adjustY(elkSection.startY))
                jsonSection.addJsonObj("startPoint", startPoint)
                
                // End Point
                val endPoint = newJsonObject
                endPoint.addJsonObj("x", edge.adjustX(elkSection.endX))
                endPoint.addJsonObj("y", edge.adjustY(elkSection.endY))
                jsonSection.addJsonObj("endPoint", endPoint)
                
                // Bend Points
                if (!elkSection.bendPoints.nullOrEmpty) {
                    val bendPoints = newJsonArray
                    elkSection.bendPoints.forEach [ pnt |
                        val jsonPnt = newJsonObject
                        jsonPnt.addJsonObj("x", edge.adjustX(pnt.x))
                        jsonPnt.addJsonObj("y", edge.adjustY(pnt.y))
                        bendPoints.addJsonArr(jsonPnt)
                    ]
                    jsonSection.addJsonObj("bendPoints", bendPoints)
                }
                
                // Incoming shape
                if (elkSection.incomingShape !== null) {
                    jsonSection.addJsonObj("incomingShape", idByElement(elkSection.incomingShape))
                }
                
                // Outgoing shape
                if (elkSection.outgoingShape !== null) {
                    jsonSection.addJsonObj("outgoingShape", idByElement(elkSection.outgoingShape))
                }
                
                // Incoming sections
                if (!elkSection.incomingSections.empty) {
                    val incomingSections = newJsonArray
                    elkSection.incomingSections.forEach [ sec |
                        incomingSections.addJsonArr(idByElement(sec))
                    ]
                    jsonSection.addJsonObj("incomingSections", incomingSections)
                }
                
                // Outgoing sections
                if (!elkSection.outgoingSections.empty) {
                    val outgoingSections = newJsonArray
                    elkSection.outgoingSections.forEach [ sec |
                        outgoingSections.addJsonArr(idByElement(sec))
                    ]
                    jsonSection.addJsonObj("outgoingSections", outgoingSections)
                }
                
            ]
            jsonObj.removeJsonObj("sections")
            jsonObj.addJsonObj("sections", sections)
        }
        
        // transfer junction points, if existent
        if (edge.hasProperty(CoreOptions.JUNCTION_POINTS)) {
            val jps = edge.getProperty(CoreOptions.JUNCTION_POINTS)
            if (!jps.nullOrEmpty) {
                val jsonJPs = newJsonArray
                jps.forEach[ jp |
                    val jsonPnt = newJsonObject
                    jsonPnt.addJsonObj("x", edge.adjustX(jp.x))
                    jsonPnt.addJsonObj("y", edge.adjustY(jp.y))
                    jsonJPs.addJsonArr(jsonPnt)
                ]
                jsonObj.addJsonObj("junctionPoints", jsonJPs)
            }
        }

        if (edge.originalParent.edgeCoordsMode === EdgeCoords.CONTAINER) {
            jsonObj.addJsonObj("container", edge.getContainingNode.identifier)
        }
    }
    
    private def dispatch transferLayoutInt2(ElkLabel label) {
        val jsonObj = labelJsonMap.get(label)
        label.recordGlobalCoords
        label.recordCoordinateModes
        label.transferShapeLayout(jsonObj)
    }
    
    private def dispatch transferLayoutInt1(Object obj) {
        // don't care about the rest
    }
    
    private def dispatch transferLayoutInt2(Object obj) {
        // don't care about the rest
    }

    private def transferShapeLayout(ElkShape shape, Object jsonObjA) {
        val jsonObj = jsonObjA.toJsonObject
        val parent = shape.jsonParent
        // pos and dimension
        jsonObj.addJsonObj("x", parent?.adjustX(shape.x) ?: shape.x)
        jsonObj.addJsonObj("y", parent?.adjustY(shape.y) ?: shape.y)
        jsonObj.addJsonObj("width", shape.width)
        jsonObj.addJsonObj("height", shape.height)
    }
    
    private def dispatch Double adjustX(ElkEdge edge, Double x) {
        val mode = edge.originalParent.edgeCoordsMode
        return switch mode {
            case EdgeCoords.ROOT: x + edge.getContainingNode.globalX
            case EdgeCoords.PARENT: x + edge.getContainingNode.globalX - edge.originalParent.globalX
            default: x
        }
    }
    
    private def dispatch Double adjustX(ElkShape shape, Double x) {
        val mode = shape.shapeCoordsMode
        return switch mode {
            case ShapeCoords.ROOT: x + shape.globalX
            default: x
        }
    }
    
    private def dispatch Double adjustY(ElkEdge edge, Double y) {
        val mode = edge.originalParent.edgeCoordsMode
        return switch mode {
            case EdgeCoords.ROOT: y + edge.getContainingNode.globalY
            case EdgeCoords.PARENT: y + edge.getContainingNode.globalY - edge.originalParent.globalY
            default: y
        }
    }
    
    private def dispatch Double adjustY(ElkShape shape, Double y) {
        val mode = shape.shapeCoordsMode
        return switch mode {
            case ShapeCoords.ROOT: y + shape.globalY
            default: y
        }
    }
    
    private def dispatch idByElement(ElkNode node) {
        return nodeIdMap.inverse.get(node)
    }
    
    private def dispatch idByElement(ElkPort port) {
        return portIdMap.inverse.get(port)
    }
    
    private def dispatch idByElement(ElkEdgeSection section) {
        return edgeSectionIdMap.inverse.get(section)
    }

    /* ---------------------------------------------------------------------------
     *                            Convenience methods
     * ---------------------------------------------------------------------------
     */
     
    private def ElkNode register(ElkNode node, Object obj) {
        val id = obj.toJsonObject.id

        nodeIdMap.put(id, node)
        nodeJsonMap.put(node, obj)

        return node
    }

    private def ElkPort register(ElkPort port, Object obj) {
        val id = obj.toJsonObject.id

        portIdMap.put(id, port)
        portJsonMap.put(port, obj)

        return port
    }

    private def ElkEdge register(ElkEdge edge, Object obj) {
        val id = obj.toJsonObject.id

        edgeIdMap.put(id, edge)
        edgeJsonMap.put(edge, obj)
        edgeOriginalParentMap.put(edge, edge.getContainingNode)

        return edge
    }

    private def ElkEdgeSection register(ElkEdgeSection edgeSection, Object obj) {
        val id = obj.toJsonObject.id

        edgeSectionIdMap.put(id, edgeSection)
        edgeSectionJsonMap.put(edgeSection, obj)

        return edgeSection
    }
    
    private def recordCoordinateModes(ElkGraphElement element) {
        val parent = element.jsonParent
        
        var scm = element.getProperty(CoreOptions.JSON_SHAPE_COORDS) ?: ShapeCoords.INHERIT
        if (scm === ShapeCoords.INHERIT) {
            scm = parent.shapeCoordsMode ?: ShapeCoords.PARENT
        }
        shapeCoordsMap.put(element, scm)
        
        var ecm = element.getProperty(CoreOptions.JSON_EDGE_COORDS) ?: EdgeCoords.INHERIT
        if (ecm === EdgeCoords.INHERIT) {
            ecm = parent.edgeCoordsMode ?: EdgeCoords.CONTAINER
        }
        edgeCoordsMap.put(element, ecm)
    }
    
    private def recordGlobalCoords(ElkShape shape) {
        val parent = shape.jsonParent
        val shapeA = parent?.shapeAncestor
        val dx = shapeA.globalX ?: 0
        val dy = shapeA.globalY ?: 0
        globalXMap.put(shape, shape.x + dx)
        globalYMap.put(shape, shape.y + dy)
    }
    
    private def dispatch shapeAncestor(ElkEdge edge) {
        return edge.originalParent
    }
    
    private def dispatch shapeAncestor(ElkShape shape) {
        return shape
    }
        
    private def dispatch jsonParent(ElkNode node) {
        return node.getParent
    }
    
    private def dispatch jsonParent(ElkPort port) {
        return port.getParent
    }
    
    private def dispatch jsonParent(ElkEdge edge) {
        return edge.originalParent
    }
    
    private def dispatch jsonParent(ElkLabel label) {
        return label.getParent
    }    
    
    private def ElkNode originalParent(ElkEdge edge) {
        return edgeOriginalParentMap.get(edge)
    }
    
    private def Double globalX(ElkShape shape) {
        return globalXMap.get(shape)
    }
    
    private def Double globalY(ElkShape shape) {
        return globalYMap.get(shape)
    }
    
    private def ShapeCoords shapeCoordsMode(ElkGraphElement element) {
        return shapeCoordsMap.get(element)
    }
    
    private def EdgeCoords edgeCoordsMode(ElkGraphElement element) {
        return edgeCoordsMap.get(element)
    }

}
