/**
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 */
package org.eclipse.elk.graph.impl;

import java.util.Map;

import org.eclipse.elk.graph.EMapPropertyHolder;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.ElkShape;

import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.IPropertyHolder;

import org.eclipse.emf.ecore.EAttribute;
import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EDataType;
import org.eclipse.emf.ecore.EGenericType;
import org.eclipse.emf.ecore.EOperation;
import org.eclipse.emf.ecore.EPackage;
import org.eclipse.emf.ecore.EReference;
import org.eclipse.emf.ecore.ETypeParameter;

import org.eclipse.emf.ecore.impl.EPackageImpl;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model <b>Package</b>.
 * <!-- end-user-doc -->
 * @generated
 */
public class ElkGraphPackageImpl extends EPackageImpl implements ElkGraphPackage {
    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass iPropertyHolderEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass eMapPropertyHolderEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkGraphElementEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkShapeEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkLabelEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkConnectableShapeEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkNodeEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkPortEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkEdgeEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkBendPointEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkEdgeSectionEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EClass elkPropertyToValueMapEntryEClass = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EDataType iPropertyEDataType = null;

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private EDataType propertyValueEDataType = null;

    /**
     * Creates an instance of the model <b>Package</b>, registered with
     * {@link org.eclipse.emf.ecore.EPackage.Registry EPackage.Registry} by the package
     * package URI value.
     * <p>Note: the correct way to create the package is via the static
     * factory method {@link #init init()}, which also performs
     * initialization of the package, or returns the registered package,
     * if one already exists.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see org.eclipse.emf.ecore.EPackage.Registry
     * @see org.eclipse.elk.graph.ElkGraphPackage#eNS_URI
     * @see #init()
     * @generated
     */
    private ElkGraphPackageImpl() {
        super(eNS_URI, ElkGraphFactory.eINSTANCE);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private static boolean isInited = false;

    /**
     * Creates, registers, and initializes the <b>Package</b> for this model, and for any others upon which it depends.
     * 
     * <p>This method is used to initialize {@link ElkGraphPackage#eINSTANCE} when that field is accessed.
     * Clients should not invoke it directly. Instead, they should simply access that field to obtain the package.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @see #eNS_URI
     * @see #createPackageContents()
     * @see #initializePackageContents()
     * @generated
     */
    public static ElkGraphPackage init() {
        if (isInited) return (ElkGraphPackage)EPackage.Registry.INSTANCE.getEPackage(ElkGraphPackage.eNS_URI);

        // Obtain or create and register package
        ElkGraphPackageImpl theElkGraphPackage = (ElkGraphPackageImpl)(EPackage.Registry.INSTANCE.get(eNS_URI) instanceof ElkGraphPackageImpl ? EPackage.Registry.INSTANCE.get(eNS_URI) : new ElkGraphPackageImpl());

        isInited = true;

        // Create package meta-data objects
        theElkGraphPackage.createPackageContents();

        // Initialize created meta-data
        theElkGraphPackage.initializePackageContents();

        // Mark meta-data to indicate it can't be changed
        theElkGraphPackage.freeze();

  
        // Update the registry and return the package
        EPackage.Registry.INSTANCE.put(ElkGraphPackage.eNS_URI, theElkGraphPackage);
        return theElkGraphPackage;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getIPropertyHolder() {
        return iPropertyHolderEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getEMapPropertyHolder() {
        return eMapPropertyHolderEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getEMapPropertyHolder_Properties() {
        return (EReference)eMapPropertyHolderEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkGraphElement() {
        return elkGraphElementEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkGraphElement_Labels() {
        return (EReference)elkGraphElementEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkGraphElement_Identifier() {
        return (EAttribute)elkGraphElementEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkShape() {
        return elkShapeEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkShape_Height() {
        return (EAttribute)elkShapeEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkShape_Width() {
        return (EAttribute)elkShapeEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkShape_X() {
        return (EAttribute)elkShapeEClass.getEStructuralFeatures().get(2);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkShape_Y() {
        return (EAttribute)elkShapeEClass.getEStructuralFeatures().get(3);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkLabel() {
        return elkLabelEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkLabel_Parent() {
        return (EReference)elkLabelEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkLabel_Text() {
        return (EAttribute)elkLabelEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkConnectableShape() {
        return elkConnectableShapeEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkConnectableShape_OutgoingEdges() {
        return (EReference)elkConnectableShapeEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkConnectableShape_IncomingEdges() {
        return (EReference)elkConnectableShapeEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkNode() {
        return elkNodeEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkNode_Ports() {
        return (EReference)elkNodeEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkNode_Children() {
        return (EReference)elkNodeEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkNode_Parent() {
        return (EReference)elkNodeEClass.getEStructuralFeatures().get(2);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkNode_ContainedEdges() {
        return (EReference)elkNodeEClass.getEStructuralFeatures().get(3);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkNode_Hierarchical() {
        return (EAttribute)elkNodeEClass.getEStructuralFeatures().get(4);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkPort() {
        return elkPortEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkPort_Parent() {
        return (EReference)elkPortEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkEdge() {
        return elkEdgeEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdge_ContainingNode() {
        return (EReference)elkEdgeEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdge_Sources() {
        return (EReference)elkEdgeEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdge_Targets() {
        return (EReference)elkEdgeEClass.getEStructuralFeatures().get(2);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdge_Sections() {
        return (EReference)elkEdgeEClass.getEStructuralFeatures().get(3);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdge_Hyperedge() {
        return (EAttribute)elkEdgeEClass.getEStructuralFeatures().get(4);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdge_Hierarchical() {
        return (EAttribute)elkEdgeEClass.getEStructuralFeatures().get(5);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdge_Selfloop() {
        return (EAttribute)elkEdgeEClass.getEStructuralFeatures().get(6);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdge_Connected() {
        return (EAttribute)elkEdgeEClass.getEStructuralFeatures().get(7);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkBendPoint() {
        return elkBendPointEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkBendPoint_X() {
        return (EAttribute)elkBendPointEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkBendPoint_Y() {
        return (EAttribute)elkBendPointEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkEdgeSection() {
        return elkEdgeSectionEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdgeSection_StartX() {
        return (EAttribute)elkEdgeSectionEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdgeSection_StartY() {
        return (EAttribute)elkEdgeSectionEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdgeSection_EndX() {
        return (EAttribute)elkEdgeSectionEClass.getEStructuralFeatures().get(2);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdgeSection_EndY() {
        return (EAttribute)elkEdgeSectionEClass.getEStructuralFeatures().get(3);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_BendPoints() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(4);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_Parent() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(5);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_OutgoingShape() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(6);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_IncomingShape() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(7);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_OutgoingSections() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(8);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EReference getElkEdgeSection_IncomingSections() {
        return (EReference)elkEdgeSectionEClass.getEStructuralFeatures().get(9);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkEdgeSection_Identifier() {
        return (EAttribute)elkEdgeSectionEClass.getEStructuralFeatures().get(10);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EClass getElkPropertyToValueMapEntry() {
        return elkPropertyToValueMapEntryEClass;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkPropertyToValueMapEntry_Key() {
        return (EAttribute)elkPropertyToValueMapEntryEClass.getEStructuralFeatures().get(0);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EAttribute getElkPropertyToValueMapEntry_Value() {
        return (EAttribute)elkPropertyToValueMapEntryEClass.getEStructuralFeatures().get(1);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EDataType getIProperty() {
        return iPropertyEDataType;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public EDataType getPropertyValue() {
        return propertyValueEDataType;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkGraphFactory getElkGraphFactory() {
        return (ElkGraphFactory)getEFactoryInstance();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private boolean isCreated = false;

    /**
     * Creates the meta-model objects for the package.  This method is
     * guarded to have no affect on any invocation but its first.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void createPackageContents() {
        if (isCreated) return;
        isCreated = true;

        // Create classes and their features
        iPropertyHolderEClass = createEClass(IPROPERTY_HOLDER);

        eMapPropertyHolderEClass = createEClass(EMAP_PROPERTY_HOLDER);
        createEReference(eMapPropertyHolderEClass, EMAP_PROPERTY_HOLDER__PROPERTIES);

        elkGraphElementEClass = createEClass(ELK_GRAPH_ELEMENT);
        createEReference(elkGraphElementEClass, ELK_GRAPH_ELEMENT__LABELS);
        createEAttribute(elkGraphElementEClass, ELK_GRAPH_ELEMENT__IDENTIFIER);

        elkShapeEClass = createEClass(ELK_SHAPE);
        createEAttribute(elkShapeEClass, ELK_SHAPE__HEIGHT);
        createEAttribute(elkShapeEClass, ELK_SHAPE__WIDTH);
        createEAttribute(elkShapeEClass, ELK_SHAPE__X);
        createEAttribute(elkShapeEClass, ELK_SHAPE__Y);

        elkLabelEClass = createEClass(ELK_LABEL);
        createEReference(elkLabelEClass, ELK_LABEL__PARENT);
        createEAttribute(elkLabelEClass, ELK_LABEL__TEXT);

        elkConnectableShapeEClass = createEClass(ELK_CONNECTABLE_SHAPE);
        createEReference(elkConnectableShapeEClass, ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES);
        createEReference(elkConnectableShapeEClass, ELK_CONNECTABLE_SHAPE__INCOMING_EDGES);

        elkNodeEClass = createEClass(ELK_NODE);
        createEReference(elkNodeEClass, ELK_NODE__PORTS);
        createEReference(elkNodeEClass, ELK_NODE__CHILDREN);
        createEReference(elkNodeEClass, ELK_NODE__PARENT);
        createEReference(elkNodeEClass, ELK_NODE__CONTAINED_EDGES);
        createEAttribute(elkNodeEClass, ELK_NODE__HIERARCHICAL);

        elkPortEClass = createEClass(ELK_PORT);
        createEReference(elkPortEClass, ELK_PORT__PARENT);

        elkEdgeEClass = createEClass(ELK_EDGE);
        createEReference(elkEdgeEClass, ELK_EDGE__CONTAINING_NODE);
        createEReference(elkEdgeEClass, ELK_EDGE__SOURCES);
        createEReference(elkEdgeEClass, ELK_EDGE__TARGETS);
        createEReference(elkEdgeEClass, ELK_EDGE__SECTIONS);
        createEAttribute(elkEdgeEClass, ELK_EDGE__HYPEREDGE);
        createEAttribute(elkEdgeEClass, ELK_EDGE__HIERARCHICAL);
        createEAttribute(elkEdgeEClass, ELK_EDGE__SELFLOOP);
        createEAttribute(elkEdgeEClass, ELK_EDGE__CONNECTED);

        elkBendPointEClass = createEClass(ELK_BEND_POINT);
        createEAttribute(elkBendPointEClass, ELK_BEND_POINT__X);
        createEAttribute(elkBendPointEClass, ELK_BEND_POINT__Y);

        elkEdgeSectionEClass = createEClass(ELK_EDGE_SECTION);
        createEAttribute(elkEdgeSectionEClass, ELK_EDGE_SECTION__START_X);
        createEAttribute(elkEdgeSectionEClass, ELK_EDGE_SECTION__START_Y);
        createEAttribute(elkEdgeSectionEClass, ELK_EDGE_SECTION__END_X);
        createEAttribute(elkEdgeSectionEClass, ELK_EDGE_SECTION__END_Y);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__BEND_POINTS);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__PARENT);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__OUTGOING_SHAPE);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__INCOMING_SHAPE);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__OUTGOING_SECTIONS);
        createEReference(elkEdgeSectionEClass, ELK_EDGE_SECTION__INCOMING_SECTIONS);
        createEAttribute(elkEdgeSectionEClass, ELK_EDGE_SECTION__IDENTIFIER);

        elkPropertyToValueMapEntryEClass = createEClass(ELK_PROPERTY_TO_VALUE_MAP_ENTRY);
        createEAttribute(elkPropertyToValueMapEntryEClass, ELK_PROPERTY_TO_VALUE_MAP_ENTRY__KEY);
        createEAttribute(elkPropertyToValueMapEntryEClass, ELK_PROPERTY_TO_VALUE_MAP_ENTRY__VALUE);

        // Create data types
        iPropertyEDataType = createEDataType(IPROPERTY);
        propertyValueEDataType = createEDataType(PROPERTY_VALUE);
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    private boolean isInitialized = false;

    /**
     * Complete the initialization of the package and its meta-model.  This
     * method is guarded to have no affect on any invocation but its first.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public void initializePackageContents() {
        if (isInitialized) return;
        isInitialized = true;

        // Initialize package
        setName(eNAME);
        setNsPrefix(eNS_PREFIX);
        setNsURI(eNS_URI);

        // Create type parameters
        addETypeParameter(iPropertyEDataType, "T");

        // Set bounds for type parameters

        // Add supertypes to classes
        eMapPropertyHolderEClass.getESuperTypes().add(this.getIPropertyHolder());
        elkGraphElementEClass.getESuperTypes().add(this.getEMapPropertyHolder());
        elkShapeEClass.getESuperTypes().add(this.getElkGraphElement());
        elkLabelEClass.getESuperTypes().add(this.getElkShape());
        elkConnectableShapeEClass.getESuperTypes().add(this.getElkShape());
        elkNodeEClass.getESuperTypes().add(this.getElkConnectableShape());
        elkPortEClass.getESuperTypes().add(this.getElkConnectableShape());
        elkEdgeEClass.getESuperTypes().add(this.getElkGraphElement());
        elkEdgeSectionEClass.getESuperTypes().add(this.getEMapPropertyHolder());

        // Initialize classes and features; add operations and parameters
        initEClass(iPropertyHolderEClass, IPropertyHolder.class, "IPropertyHolder", IS_ABSTRACT, IS_INTERFACE, !IS_GENERATED_INSTANCE_CLASS);

        EOperation op = addEOperation(iPropertyHolderEClass, this.getIPropertyHolder(), "setProperty", 0, 1, IS_UNIQUE, IS_ORDERED);
        ETypeParameter t1 = addETypeParameter(op, "T");
        EGenericType g1 = createEGenericType(this.getIProperty());
        EGenericType g2 = createEGenericType();
        g1.getETypeArguments().add(g2);
        EGenericType g3 = createEGenericType(t1);
        g2.setELowerBound(g3);
        addEParameter(op, g1, "property", 0, 1, IS_UNIQUE, IS_ORDERED);
        g1 = createEGenericType(t1);
        addEParameter(op, g1, "value", 0, 1, IS_UNIQUE, IS_ORDERED);

        op = addEOperation(iPropertyHolderEClass, null, "getProperty", 0, 1, IS_UNIQUE, IS_ORDERED);
        t1 = addETypeParameter(op, "T");
        g1 = createEGenericType(this.getIProperty());
        g2 = createEGenericType(t1);
        g1.getETypeArguments().add(g2);
        addEParameter(op, g1, "property", 0, 1, IS_UNIQUE, IS_ORDERED);
        g1 = createEGenericType(t1);
        initEOperation(op, g1);

        op = addEOperation(iPropertyHolderEClass, ecorePackage.getEBoolean(), "hasProperty", 0, 1, IS_UNIQUE, IS_ORDERED);
        g1 = createEGenericType(this.getIProperty());
        g2 = createEGenericType();
        g1.getETypeArguments().add(g2);
        addEParameter(op, g1, "property", 0, 1, IS_UNIQUE, IS_ORDERED);

        op = addEOperation(iPropertyHolderEClass, this.getIPropertyHolder(), "copyProperties", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, this.getIPropertyHolder(), "source", 0, 1, IS_UNIQUE, IS_ORDERED);

        op = addEOperation(iPropertyHolderEClass, null, "getAllProperties", 0, 1, IS_UNIQUE, IS_ORDERED);
        g1 = createEGenericType(ecorePackage.getEMap());
        g2 = createEGenericType(this.getIProperty());
        g1.getETypeArguments().add(g2);
        g3 = createEGenericType();
        g2.getETypeArguments().add(g3);
        g2 = createEGenericType(ecorePackage.getEJavaObject());
        g1.getETypeArguments().add(g2);
        initEOperation(op, g1);

        initEClass(eMapPropertyHolderEClass, EMapPropertyHolder.class, "EMapPropertyHolder", IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getEMapPropertyHolder_Properties(), this.getElkPropertyToValueMapEntry(), null, "properties", null, 0, -1, EMapPropertyHolder.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        initEClass(elkGraphElementEClass, ElkGraphElement.class, "ElkGraphElement", IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkGraphElement_Labels(), this.getElkLabel(), this.getElkLabel_Parent(), "labels", null, 0, -1, ElkGraphElement.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkGraphElement_Identifier(), ecorePackage.getEString(), "identifier", null, 0, 1, ElkGraphElement.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        initEClass(elkShapeEClass, ElkShape.class, "ElkShape", IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEAttribute(getElkShape_Height(), ecorePackage.getEDouble(), "height", "0.0", 1, 1, ElkShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkShape_Width(), ecorePackage.getEDouble(), "width", "0.0", 1, 1, ElkShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkShape_X(), ecorePackage.getEDouble(), "x", "0.0", 1, 1, ElkShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkShape_Y(), ecorePackage.getEDouble(), "y", "0.0", 1, 1, ElkShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        op = addEOperation(elkShapeEClass, null, "setDimensions", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "width", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "height", 0, 1, IS_UNIQUE, IS_ORDERED);

        op = addEOperation(elkShapeEClass, null, "setLocation", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "x", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "y", 0, 1, IS_UNIQUE, IS_ORDERED);

        initEClass(elkLabelEClass, ElkLabel.class, "ElkLabel", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkLabel_Parent(), this.getElkGraphElement(), this.getElkGraphElement_Labels(), "parent", null, 0, 1, ElkLabel.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkLabel_Text(), ecorePackage.getEString(), "text", "", 0, 1, ElkLabel.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        initEClass(elkConnectableShapeEClass, ElkConnectableShape.class, "ElkConnectableShape", IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkConnectableShape_OutgoingEdges(), this.getElkEdge(), this.getElkEdge_Sources(), "outgoingEdges", null, 0, -1, ElkConnectableShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkConnectableShape_IncomingEdges(), this.getElkEdge(), this.getElkEdge_Targets(), "incomingEdges", null, 0, -1, ElkConnectableShape.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        initEClass(elkNodeEClass, ElkNode.class, "ElkNode", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkNode_Ports(), this.getElkPort(), this.getElkPort_Parent(), "ports", null, 0, -1, ElkNode.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkNode_Children(), this.getElkNode(), this.getElkNode_Parent(), "children", null, 0, -1, ElkNode.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkNode_Parent(), this.getElkNode(), this.getElkNode_Children(), "parent", null, 0, 1, ElkNode.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkNode_ContainedEdges(), this.getElkEdge(), this.getElkEdge_ContainingNode(), "containedEdges", null, 0, -1, ElkNode.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkNode_Hierarchical(), ecorePackage.getEBoolean(), "hierarchical", null, 0, 1, ElkNode.class, IS_TRANSIENT, IS_VOLATILE, !IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, IS_DERIVED, IS_ORDERED);

        initEClass(elkPortEClass, ElkPort.class, "ElkPort", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkPort_Parent(), this.getElkNode(), this.getElkNode_Ports(), "parent", null, 0, 1, ElkPort.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        initEClass(elkEdgeEClass, ElkEdge.class, "ElkEdge", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEReference(getElkEdge_ContainingNode(), this.getElkNode(), this.getElkNode_ContainedEdges(), "containingNode", null, 0, 1, ElkEdge.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdge_Sources(), this.getElkConnectableShape(), this.getElkConnectableShape_OutgoingEdges(), "sources", null, 0, -1, ElkEdge.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdge_Targets(), this.getElkConnectableShape(), this.getElkConnectableShape_IncomingEdges(), "targets", null, 0, -1, ElkEdge.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdge_Sections(), this.getElkEdgeSection(), this.getElkEdgeSection_Parent(), "sections", null, 0, -1, ElkEdge.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdge_Hyperedge(), ecorePackage.getEBoolean(), "hyperedge", null, 0, 1, ElkEdge.class, IS_TRANSIENT, IS_VOLATILE, !IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdge_Hierarchical(), ecorePackage.getEBoolean(), "hierarchical", null, 0, 1, ElkEdge.class, IS_TRANSIENT, IS_VOLATILE, !IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdge_Selfloop(), ecorePackage.getEBoolean(), "selfloop", null, 0, 1, ElkEdge.class, IS_TRANSIENT, IS_VOLATILE, !IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdge_Connected(), ecorePackage.getEBoolean(), "connected", null, 0, 1, ElkEdge.class, IS_TRANSIENT, IS_VOLATILE, !IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, IS_DERIVED, IS_ORDERED);

        initEClass(elkBendPointEClass, ElkBendPoint.class, "ElkBendPoint", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEAttribute(getElkBendPoint_X(), ecorePackage.getEDouble(), "x", "0.0", 1, 1, ElkBendPoint.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkBendPoint_Y(), ecorePackage.getEDouble(), "y", "0.0", 1, 1, ElkBendPoint.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        op = addEOperation(elkBendPointEClass, null, "set", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "x", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "y", 0, 1, IS_UNIQUE, IS_ORDERED);

        initEClass(elkEdgeSectionEClass, ElkEdgeSection.class, "ElkEdgeSection", !IS_ABSTRACT, !IS_INTERFACE, IS_GENERATED_INSTANCE_CLASS);
        initEAttribute(getElkEdgeSection_StartX(), ecorePackage.getEDouble(), "startX", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdgeSection_StartY(), ecorePackage.getEDouble(), "startY", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdgeSection_EndX(), ecorePackage.getEDouble(), "endX", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdgeSection_EndY(), ecorePackage.getEDouble(), "endY", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_BendPoints(), this.getElkBendPoint(), null, "bendPoints", null, 0, -1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_Parent(), this.getElkEdge(), this.getElkEdge_Sections(), "parent", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, !IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_OutgoingShape(), this.getElkConnectableShape(), null, "outgoingShape", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_IncomingShape(), this.getElkConnectableShape(), null, "incomingShape", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_OutgoingSections(), this.getElkEdgeSection(), this.getElkEdgeSection_IncomingSections(), "outgoingSections", null, 0, -1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEReference(getElkEdgeSection_IncomingSections(), this.getElkEdgeSection(), this.getElkEdgeSection_OutgoingSections(), "incomingSections", null, 0, -1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_COMPOSITE, IS_RESOLVE_PROXIES, !IS_UNSETTABLE, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkEdgeSection_Identifier(), ecorePackage.getEString(), "identifier", null, 0, 1, ElkEdgeSection.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        op = addEOperation(elkEdgeSectionEClass, null, "setStartLocation", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "x", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "y", 0, 1, IS_UNIQUE, IS_ORDERED);

        op = addEOperation(elkEdgeSectionEClass, null, "setEndLocation", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "x", 0, 1, IS_UNIQUE, IS_ORDERED);
        addEParameter(op, ecorePackage.getEDouble(), "y", 0, 1, IS_UNIQUE, IS_ORDERED);

        initEClass(elkPropertyToValueMapEntryEClass, Map.Entry.class, "ElkPropertyToValueMapEntry", !IS_ABSTRACT, !IS_INTERFACE, !IS_GENERATED_INSTANCE_CLASS);
        g1 = createEGenericType(this.getIProperty());
        g2 = createEGenericType();
        g1.getETypeArguments().add(g2);
        initEAttribute(getElkPropertyToValueMapEntry_Key(), g1, "key", null, 0, 1, Map.Entry.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);
        initEAttribute(getElkPropertyToValueMapEntry_Value(), this.getPropertyValue(), "value", null, 0, 1, Map.Entry.class, !IS_TRANSIENT, !IS_VOLATILE, IS_CHANGEABLE, !IS_UNSETTABLE, !IS_ID, IS_UNIQUE, !IS_DERIVED, IS_ORDERED);

        // Initialize data types
        initEDataType(iPropertyEDataType, IProperty.class, "IProperty", IS_SERIALIZABLE, !IS_GENERATED_INSTANCE_CLASS);
        initEDataType(propertyValueEDataType, Object.class, "PropertyValue", IS_SERIALIZABLE, !IS_GENERATED_INSTANCE_CLASS);

        // Create resource
        createResource(eNS_URI);
    }

} //ElkGraphPackageImpl
