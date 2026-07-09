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

import org.eclipse.elk.graph.*;

import org.eclipse.elk.graph.properties.IProperty;

import org.eclipse.emf.ecore.EClass;
import org.eclipse.emf.ecore.EDataType;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EPackage;

import org.eclipse.emf.ecore.impl.EFactoryImpl;

import org.eclipse.emf.ecore.plugin.EcorePlugin;

/**
 * <!-- begin-user-doc -->
 * An implementation of the model <b>Factory</b>.
 * <!-- end-user-doc -->
 * @generated
 */
public class ElkGraphFactoryImpl extends EFactoryImpl implements ElkGraphFactory {
    /**
     * Creates the default factory implementation.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public static ElkGraphFactory init() {
        try {
            ElkGraphFactory theElkGraphFactory = (ElkGraphFactory)EPackage.Registry.INSTANCE.getEFactory(ElkGraphPackage.eNS_URI);
            if (theElkGraphFactory != null) {
                return theElkGraphFactory;
            }
        }
        catch (Exception exception) {
            EcorePlugin.INSTANCE.log(exception);
        }
        return new ElkGraphFactoryImpl();
    }

    /**
     * Creates an instance of the factory.
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkGraphFactoryImpl() {
        super();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public EObject create(EClass eClass) {
        switch (eClass.getClassifierID()) {
            case ElkGraphPackage.ELK_LABEL: return createElkLabel();
            case ElkGraphPackage.ELK_NODE: return createElkNode();
            case ElkGraphPackage.ELK_PORT: return createElkPort();
            case ElkGraphPackage.ELK_EDGE: return createElkEdge();
            case ElkGraphPackage.ELK_BEND_POINT: return createElkBendPoint();
            case ElkGraphPackage.ELK_EDGE_SECTION: return createElkEdgeSection();
            case ElkGraphPackage.ELK_PROPERTY_TO_VALUE_MAP_ENTRY: return (EObject)createElkPropertyToValueMapEntry();
            default:
                throw new IllegalArgumentException("The class '" + eClass.getName() + "' is not a valid classifier");
        }
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public Object createFromString(EDataType eDataType, String initialValue) {
        switch (eDataType.getClassifierID()) {
            case ElkGraphPackage.IPROPERTY:
                return createIPropertyFromString(eDataType, initialValue);
            case ElkGraphPackage.PROPERTY_VALUE:
                return createPropertyValueFromString(eDataType, initialValue);
            default:
                throw new IllegalArgumentException("The datatype '" + eDataType.getName() + "' is not a valid classifier");
        }
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    @Override
    public String convertToString(EDataType eDataType, Object instanceValue) {
        switch (eDataType.getClassifierID()) {
            case ElkGraphPackage.IPROPERTY:
                return convertIPropertyToString(eDataType, instanceValue);
            case ElkGraphPackage.PROPERTY_VALUE:
                return convertPropertyValueToString(eDataType, instanceValue);
            default:
                throw new IllegalArgumentException("The datatype '" + eDataType.getName() + "' is not a valid classifier");
        }
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkLabel createElkLabel() {
        ElkLabelImpl elkLabel = new ElkLabelImpl();
        return elkLabel;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkNode createElkNode() {
        ElkNodeImpl elkNode = new ElkNodeImpl();
        return elkNode;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkPort createElkPort() {
        ElkPortImpl elkPort = new ElkPortImpl();
        return elkPort;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkEdge createElkEdge() {
        ElkEdgeImpl elkEdge = new ElkEdgeImpl();
        return elkEdge;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkBendPoint createElkBendPoint() {
        ElkBendPointImpl elkBendPoint = new ElkBendPointImpl();
        return elkBendPoint;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkEdgeSection createElkEdgeSection() {
        ElkEdgeSectionImpl elkEdgeSection = new ElkEdgeSectionImpl();
        return elkEdgeSection;
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public Map.Entry<IProperty<?>, Object> createElkPropertyToValueMapEntry() {
        ElkPropertyToValueMapEntryImpl elkPropertyToValueMapEntry = new ElkPropertyToValueMapEntryImpl();
        return elkPropertyToValueMapEntry;
    }

    /**
     * <!-- begin-user-doc -->
     * Deserialization requires the layout meta data service to be present and can thus not happen in this plug-in.
     * This method thus always returns {@code null}.
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public IProperty<?> createIPropertyFromString(EDataType eDataType, String initialValue) {
        return null;
    }

    /**
     * <!-- begin-user-doc -->
     * Serializes {@link IProperty} instances to their ID. This is used while serializing {@link EMapPropertyHolder}s.
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public String convertIPropertyToString(EDataType eDataType, Object instanceValue) {
        return ((IProperty<?>) instanceValue).getId();
    }

    /**
     * <!-- begin-user-doc -->
     * Deserialization requires the layout meta data service to be present and can thus not happen in this plug-in.
     * This method thus always returns {@code null}.
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public Object createPropertyValueFromString(EDataType eDataType, String initialValue) {
        return null;
    }

    /**
     * <!-- begin-user-doc -->
     * Serializes property values to their String representations by calling {@link Object#toString()}. This is used
     * while serializing {@link EMapPropertyHolder}s. If a property value implements {@code IDataObject} it will be
     * able to parse this String representation during deserialization.
     * <!-- end-user-doc -->
     * @generated NOT
     */
    public String convertPropertyValueToString(EDataType eDataType, Object instanceValue) {
        return instanceValue.toString();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @generated
     */
    public ElkGraphPackage getElkGraphPackage() {
        return (ElkGraphPackage)getEPackage();
    }

    /**
     * <!-- begin-user-doc -->
     * <!-- end-user-doc -->
     * @deprecated
     * @generated
     */
    @Deprecated
    public static ElkGraphPackage getPackage() {
        return ElkGraphPackage.eINSTANCE;
    }

} //ElkGraphFactoryImpl
